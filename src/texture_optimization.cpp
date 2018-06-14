#include "texture_optimization.h"

#include "uv.h"
#include "mesh_graph.h"
#include "timer.h"
#include "dcp_solver.h"
#include "uniform_solver.h"
#include "iterative.h"
#include "metric.h"
#include "parameterization_checker.h"
#include "mesh_utils.h"

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/complex/algorithms/crease_cut.h>
#include <vcg/space/rasterized_outline2_packer.h>
#include <vcg/space/outline2_packer.h>
#include <wrap/qt/outline2_rasterizer.h>
#include <wrap/qt/Outline2ToQImage.h>
#include <vcg/space/index/grid_util2d.h>
#include <vcg/space/segment2.h>
#include <vcg/space/intersection2.h>

#include <wrap/io_trimesh/export.h>

//#include<vcg/complex/algorithms/cut_tree.h>
//#include<vcg/complex/algorithms/curve_on_manifold.h>
//#include<vcg/complex/algorithms/crease_cut.h>

#include <vector>
#include <algorithm>

static bool SegmentBoxIntersection(const Segment2<double>& seg, const Box2d& box)
{
    Point2d isec;
    Point2d c1{box.min};
    Point2d c2{box.max[0], box.min[1]};
    Point2d c3{box.max};
    Point2d c4{box.min[0], box.max[1]};

    if (SegmentSegmentIntersection(seg, Segment2<double>{c1, c2}, isec)) return true;
    if (SegmentSegmentIntersection(seg, Segment2<double>{c2, c3}, isec)) return true;
    if (SegmentSegmentIntersection(seg, Segment2<double>{c3, c4}, isec)) return true;
    if (SegmentSegmentIntersection(seg, Segment2<double>{c4, c2}, isec)) return true;
    if (SegmentSegmentIntersection(seg, Segment2<double>{c1, c2}, isec)) return true;

    // if the segment does not intersect the sides, check if it is fully contained in the box
    return (box.min[0] <= std::min(seg.P0()[0], seg.P1()[0]) &&
            box.min[1] <= std::min(seg.P0()[1], seg.P1()[1]) &&
            box.max[0] >= std::max(seg.P0()[0], seg.P1()[0]) &&
            box.max[1] >= std::max(seg.P0()[1], seg.P1()[1]));
}

bool ChartParameterizationHasOverlaps(Mesh& m, GraphManager::ChartHandle chart)
{
    using OutlineIndex = std::pair<std::size_t, std::size_t>;

    std::vector<std::vector<Point2d>> outlines;
    ChartOutlinesUV(m, chart, outlines);

    std::unordered_map<Point2i, std::vector<OutlineIndex>, Point2iHasher> grid;

    std::size_t elems = 0;
    for (auto& o : outlines) {
        elems += o.size();
    }

    // init grid helper
    BasicGrid2D<double> gh;
    gh.bbox = chart->UVBox();
    BestDim2D<double>(elems, gh.bbox.Dim(), gh.siz);
    gh.ComputeDimAndVoxel();

    // populate grid with edges
    for (std::size_t i = 0; i < outlines.size(); ++i) {
        auto sz = outlines[i].size();
        for (std::size_t j = 0; j < sz; ++j) {
            Point2d a = outlines[i][j];
            Point2d b = outlines[i][(j+1)%sz];
            Box2d edgeBox;
            edgeBox.Add(a);
            edgeBox.Add(b);
            Box2i gridCover;
            gh.BoxToIBox(edgeBox, gridCover);
            Segment2<double> e{a, b};
            for (int h = gridCover.min[0]; h <= gridCover.max[0]; h++) {
                for (int k = gridCover.min[1]; k <= gridCover.max[1]; k++) {
                    Box2d cell;
                    Point2i voxel{h, k};
                    gh.IPiToBox(voxel, cell);
                    if (SegmentBoxIntersection(e, cell)) {
                        grid[voxel].push_back(std::make_pair(i, j));
                    }
                }
            }
        }
    }

    for (auto& entry : grid) {
        for (auto i1 : entry.second) {
            for (auto i2 : entry.second) {
                if (i1 == i2) continue;
                auto sz1 = outlines[i1.first].size();
                Point2d a1 = outlines[i1.first][i1.second];
                Point2d b1 = outlines[i1.first][(i1.second + 1) % sz1];

                auto sz2 = outlines[i2.first].size();
                Point2d a2 = outlines[i2.first][i2.second];
                Point2d b2 = outlines[i2.first][(i2.second + 1) % sz2];

                Segment2<double> s1{a1, b1};
                Segment2<double> s2{a2, b2};
                Point2d intersectionPoint;
                if (SegmentSegmentIntersection(s1, s2, intersectionPoint)
                        && intersectionPoint != a1 && intersectionPoint != b1
                        && intersectionPoint != a2 && intersectionPoint != b2) {
                    std::cout << "( " <<  a1[0] << " , " << a1[1] << " ) (" << b1[0] << " , " << b1[1] << " )" << std::endl;
                    std::cout << "( " <<  a2[0] << " , " << a2[1] << " ) (" << b2[0] << " , " << b2[1] << " )" << std::endl;
                    std::cout << intersectionPoint[0] << " , " << intersectionPoint[1] << std::endl;
                    return true;
                }
            }
        }
    }

    return false;
}

void PreprocessMesh(Mesh& m)
{
    std::cout << "FIXME" << std::endl;
    // Compute preliminary parameterization graph
    auto graph = ComputeParameterizationGraph(m, nullptr);
    // Parameterize regions that are not parameterized
    ReparameterizeZeroAreaRegions(m, graph);
}

void ReparameterizeZeroAreaRegions(Mesh &m, std::shared_ptr<MeshGraph> graph)
{
    assert(HasParameterizationScaleInfoAttribute(m));
    auto info = GetParameterizationScaleInfoAttribute(m);
    double scale = info().scale;

    int numNoParam = 0;
    int numParameterized = 0;

    ParameterizationStrategy strategy;
    strategy.directParameterizer = FixedBorderBijective;
    strategy.energy = EnergyType::SymmetricDirichlet;
    strategy.geometry = Model;
    strategy.descent = ScalableLocallyInjectiveMappings;
    strategy.optimizerIterations = 200;
    strategy.padBoundaries = true;
    strategy.applyCut = false;
    strategy.warmStart = false;

    for (auto& entry : graph->charts) {
        auto chart = entry.second;
        if (chart->AreaUV() > 0)
            continue;
        numNoParam++;

        std::cout << "Parameterizing region of " << chart->FN() << " zero UV area faces" << std::endl;

        ParameterizerObject po{chart, strategy};
        bool parameterized = po.Parameterize();

        if (!parameterized) {
            std::cout << "WARNING: preliminary parameterization of chart " << chart->id << " failed" << std::endl;
        } else {
            po.SyncChart();
            // As convenience to detect regions that originally did not have a parameterization, the texcoords
            // of such regions have negative u and a randomly displaced v
            Box2d box = chart->UVBox();
            double randomDisplacementV = rand() / (double) RAND_MAX;
            for (auto fptr : chart->fpVec) {
                for (int i = 0; i < 3; ++i) {
                    fptr->WT(i).P().X() -= (box.min.X() + box.DimX());
                    fptr->WT(i).P().Y() -= (box.min.Y() + randomDisplacementV);
                    fptr->WT(i).P() *= scale;
                }
            }
            numParameterized++;
        }
    }

    std::cout << "[LOG] Newly parameterized regions: " << numParameterized << "/" << numNoParam << std::endl;
}

/* Ugly function that tries to perform subsequent merges after a split. It uses
 * a very simple heuristic that selects the two largest charts in the split, and
 * iteratively grows them until all the charts in the split are covered,
 * producing to 2 new aggregates that will be parameterized independently */
static void RecoverFromSplit(std::vector<ChartHandle>& split, GraphManager& gm, std::vector<ChartHandle>& chartVec)
{
    chartVec.clear();

    std::sort(split.begin(), split.end(),
              [](const ChartHandle& c1, const ChartHandle& c2) { return c1->Area3D() > c2->Area3D(); }
    );
    ChartHandle c1 = split[0];
    ChartHandle c2 = split[1];

    std::unordered_set<ChartHandle, FaceGroup::Hasher> charts;
    for (std::size_t i = 2; i < split.size(); ++i) {
        charts.insert(split[i]);
    }

    while (charts.size() > 0) {
        // grow c1
        for (auto ch : c1->adj) {
            if (charts.count(ch) == 1) {
                charts.erase(ch);
                c1 = gm.Collapse(GraphManager::Edge(c1, ch));
                break;
            }
        }

        // grow c2
        for (auto ch : c2->adj) {
            if (charts.count(ch) == 1) {
                charts.erase(ch);
                c2 = gm.Collapse(GraphManager::Edge(c2, ch));
                break;
            }
        }
    }

    chartVec.push_back(c1);
    chartVec.push_back(c2);
    std::cout << "Recovery produced two charts of sizes " << c1->numMerges + 1 << " (" << c1->id << ") and "
              << c2->numMerges + 1 << " (" << c2->id << ")" << std::endl;
}

int ParameterizeGraph(GraphManager& gm, ParameterizationStrategy strategy, double injectivityTolerance, bool retry)
{
    struct ParamTask {
        ChartHandle chart;
        bool warmStart;
    };

    bool injectivityCheckRequired;
    if (injectivityTolerance < 0)
        injectivityCheckRequired = false;
    else {
        assert(injectivityTolerance >= 0 && injectivityTolerance <= 1);
        injectivityCheckRequired = true;
    }

    Timer timer;

    GraphHandle graph = gm.Graph();

    assert(HasWedgeTexCoordStorageAttribute(graph->mesh));
    auto wtcsattr = GetWedgeTexCoordStorageAttribute(graph->mesh);

    std::cout << "Parameterizing " << graph->charts.size() << " regions..." << std::endl;

    // gather the charts to parameterize
    std::deque<ParamTask> paramQueue;
    for (auto entry : graph->charts) {
        ChartHandle chart = entry.second;
        if (chart->numMerges > 0)
            paramQueue.push_back(ParamTask{chart, false});
    }

    int iter = 0;
    int numFailed = 0;
    while (paramQueue.size() > 0) {
        ParamTask task = paramQueue.front();
        paramQueue.pop_front();

        ChartHandle chart = task.chart;
        std::cout << "Chart " << chart->id << " - FN=" << chart->FN() << ", FI=" << tri::Index(graph->mesh, chart->Fp()) << std::endl;
        if (chart->numMerges == 0) {
            // if there are no merges and the ws flag is set, then this chart
            // is the result of a final split, so restore the initial texture coordinates
            if (task.warmStart == true) {
                for (auto fptr : chart->fpVec) {
                    TexCoordStorage tcs = wtcsattr[fptr];
                    for (int i = 0; i < 3; ++i) {
                        fptr->WT(i) = tcs.tc[i];
                        fptr->V(i)->T() = tcs.tc[i];
                    }
                }
            }
            continue;
        }

        ParameterizationStrategy strat = strategy;
        if (task.warmStart)
            strat.warmStart = true;

        ParameterizerObject po{chart, strat};
        bool ok = po.Parameterize();
        if (ok) {
            std::cout << "  Chart parameterized correctly" << std::endl;
            po.SyncChart();
            if (injectivityCheckRequired) {
                RasterizedParameterizationStats stats = GetRasterizationStats(chart, 1024, 1024);
                double fraction = stats.lostFragments / (double) stats.totalFragments;
                if (fraction > injectivityTolerance) {
                    std::cout << "WARNING: REGION " << chart->id << " HAS OVERLAPS IN THE PARAMETERIZATION "
                              << "(overlap fraction = " << fraction << ")" << std::endl;
                    ok = false;
                }
            }
        } else {
            numFailed++;
            std::cout << "Parameterization failed" << std::endl;
        }

        if (ok) {
            // If the parameterization is valid, scale chart relative to the original parameterization area value
            // Normalize area: the region gets scaled by sqrt(oldUVArea)/sqrt(newUVArea) to keep the original proportions
            // between the regions of the atlas in an attempt to not lose too much texture data when the new texture is rendered
            double oldUvArea = chart->OriginalAreaUV();
            double newUvArea = chart->AreaUV();
            double scale = std::sqrt(oldUvArea / newUvArea);
            assert(scale > 0);
            // scale shoud be very close to 1 if we optimize for area distortion wrt to original uv coords
            std::cout << "Chart scale value = " << scale << std::endl;
            vcg::Box2d uvBox = chart->UVBox();
            for (auto fptr : chart->fpVec) {
                for (int i = 0; i < 3; ++i) {
                    fptr->WT(i).P() = (fptr->WT(i).P() - uvBox.min) * scale;
                }
            }
            chart->numMerges = 0;
        } else {
            // split the chart in the mesh graph using the graph manager, and readd
            // the new task with the warmStart flag set to true
            std::vector<ChartHandle> splitCharts;
            gm.Split(chart->id, splitCharts);
            if (retry) {
                std::vector<ChartHandle> newCharts;
                RecoverFromSplit(splitCharts, gm, newCharts);
                for (auto& c : newCharts)
                    paramQueue.push_back(ParamTask{c, true});
                    //paramQueue.push_back(ParamTask{c, false});
            } else {
                // restore original texture coordinates for each chart
                for (auto split : splitCharts) {
                    for (auto fptr : split->fpVec) {
                        TexCoordStorage tcs = wtcsattr[fptr];
                        for (int i = 0; i < 3; ++i) {
                            fptr->WT(i) = tcs.tc[i];
                            fptr->V(i)->T() = tcs.tc[i];
                        }
                    }
                }
            }
        }

        std::cout << "Iteration " << iter++ << " took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    }

    // Pack the atlas

    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(graph->mesh);

    std::vector<std::vector<Point2f>> texOutlines;
    std::unordered_map<RegionID,std::size_t> outlineMap; // map each region to the index of its outline in texOutlines

    for (auto entry : graph->charts) {
        GraphManager::ChartHandle chart = entry.second;
        std::cout << "Chart " << chart->id << " - FN=" << chart->FN() << ", FI=" << tri::Index(graph->mesh, chart->Fp()) << std::endl;

        chart->ParameterizationChanged(); // packing changes texture coords even for charts that are not reparameterized

        // Save the outline of the parameterization for this portion of the mesh
        std::vector<std::vector<Point2f>> uvOutlines;
        ChartOutlinesUV(graph->mesh, chart, uvOutlines);
        int i = tri::OutlineUtil<float>::LargestOutline2(uvOutlines);
        if (tri::OutlineUtil<float>::Outline2Area(uvOutlines[i]) < 0)
            tri::OutlineUtil<float>::ReverseOutline2(uvOutlines[i]);

        outlineMap[chart->id] = texOutlines.size();
        texOutlines.push_back(uvOutlines[i]);
    }

    std::cout << "Parameterization took " << timer.TimeElapsed() << " seconds" << std::endl;

    std::cout << "Packing the atlas..." << std::endl;

    // pack the atlas TODO function parameter to choose the packing strategy
    RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters packingParam;
    packingParam.costFunction  = RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters::LowestHorizon;
    packingParam.doubleHorizon = true;
    packingParam.cellSize = 2;
    packingParam.pad = 4;
    packingParam.rotationNum = 16; //number of rasterizations in 90°

    TextureObjectHandle to = gm.Graph()->textureObject;
    int gridWidth = to->TextureWidth(0) / 2;
    int gridHeight = to->TextureHeight(0) / 2;
    Point2i gridSize(gridWidth, gridHeight);
    std::vector<Similarity2f> transforms;

    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(graph->mesh);

    RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Pack(texOutlines, gridSize, transforms, packingParam);

    Point2f cover;
    //PolyPacker<float>::PackAsObjectOrientedRect(texOutlines, gridSize, transforms, cover);
    //PolyPacker<float>::PackAsAxisAlignedRect(texOutlines, gridSize, transforms, cover);

    std::cout << "Packing took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;

    //assert(transforms.size() == pdata.charts.size());

    for (auto p : outlineMap) {
        for (auto fptr : graph->charts[p.first]->fpVec) {
            for (int j = 0; j < fptr->VN(); ++j) {
                Point2d uv = fptr->WT(j).P();
                Point2f transformedTexCoordPos = transforms[p.second] * (Point2f(uv[0], uv[1]));
                fptr->WT(j).P() = Point2d{transformedTexCoordPos[0] / double(gridSize.X()), transformedTexCoordPos[1] / double(gridSize.Y())};
            }
        }
    }

    return numFailed;
}

void RecomputeSegmentation(GraphManager &gm, std::size_t regionCount, std::size_t minRegionSize)
{
    Timer timer;
    int mergeCount;
    int numIter = 0;

    assert(regionCount > 0);

    std::cout << "[LOG] Reduction strategy TargetRegionCount=" << regionCount << " (region threshold " << minRegionSize << ")" << std::endl;

    do {
        mergeCount = gm.CloseMacroRegions(minRegionSize);
        //mergeCount = 0;

        while (gm.HasNextEdge()) {
            auto we = gm.PeekNextEdge();
            bool regionReached = (regionCount <= 0) || (gm.Graph()->Count() <= regionCount);
            bool sizeThresholdReached = (we.first.a->FN() > minRegionSize && we.first.b->FN() > minRegionSize);
            if (regionReached && sizeThresholdReached)
                break;
            else {
                gm.RemoveNextEdge();
                gm.Collapse(we.first);
                mergeCount++;
                if (mergeCount%50 == 0) {
                    std::cout << "Merged " << mergeCount << " regions..." << std::endl;
                }
            }
        }

        std::cout << "Iteration "  << numIter << " took " << timer.TimeSinceLastCheck() << " seconds ("
                  << mergeCount << " merges took place)" << std::endl;

        numIter++;

    } while (mergeCount > 0);

    std::cout << "Stopping after " << numIter << " passes and " << timer.TimeElapsed() << " seconds" << std::endl;
}
