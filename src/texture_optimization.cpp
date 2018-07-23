#include "texture_optimization.h"

#include "mesh.h"
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

bool ChartParameterizationHasOverlaps(Mesh& m, FaceGroup& chart)
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
    gh.bbox = chart.UVBox();
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
                    //goto isec_detected;
                }
            }
        }
    }

    return false;

    isec_detected:

    static int n = 0;

    MyMesh edgeMesh;

    std::vector<std::vector<Point3d>> outline3Vec;
    for (auto& outline : outlines) {
        std::vector<Point3d> vec;
        for (auto& p : outline) {
            vec.push_back(Point3d(p.X(), p.Y(), 0));
        }
        outline3Vec.push_back(vec);
    }
    tri::OutlineUtil<double>::ConvertOutline3VecToEdgeMesh(outline3Vec, edgeMesh);

    std::stringstream ss;
    ss << "outline_intersect_" << n++ << ".obj";

    tri::io::Exporter<MyMesh>::Save(edgeMesh, ss.str().c_str(), io::Mask::IOM_VERTCOORD);

    return true;

}

void RecoverFromSplit(std::vector<ChartHandle>& split, GraphManager& gm, std::vector<ChartHandle>& chartVec, bool binarySplit)
{
    chartVec.clear();

    GraphHandle graph = gm.Graph();
    Mesh& m = graph->mesh;
    auto ccid = GetConnectedComponentIDAttribute(m);

    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(m);

    std::vector<Mesh::FacePointer> faces;
    for (auto ch : split) {
        for (auto fptr : ch->fpVec) faces.push_back(fptr);
    }

    std::vector<std::unordered_set<ChartHandle, FaceGroup::Hasher>> mergeSets;

    std::unordered_set<Mesh::FacePointer> visited;

    for (auto fptr : faces) {
        if (visited.count(fptr) == 0) {
            std::unordered_set<ChartHandle, FaceGroup::Hasher> currentMerge;
            std::unordered_set<RegionID> forbidden;
            FaceGroup aggregate{m, INVALID_ID};
            std::stack<Mesh::FacePointer> s;
            std::unordered_set<Mesh::FacePointer> inStack;
            s.push(fptr);
            inStack.insert(fptr);
            //std::cout << "Seed face " << tri::Index(m, fptr) << std::endl;
            while (!s.empty()) {
                Mesh::FacePointer fp = s.top();
                s.pop();
                assert(visited.count(fp) == 0);

                assert(inStack.count(fp) > 0);
                inStack.erase(fp);

                RegionID id = ccid[fp];
                ChartHandle chfp = graph->GetChart(id);
                bool expandVisit = false;

                if (forbidden.count(id) > 0)
                    continue;

                if (currentMerge.count(chfp) > 0) {
                    expandVisit = true;
                } else {
                    // see if adding the new region causes overlaps
                    for (auto ffp : chfp->fpVec)
                        aggregate.AddFace(ffp);

                    if (ChartParameterizationHasOverlaps(m, aggregate)) {
                        if (currentMerge.size() == 0) {
                            // handle special case where a single subchart is overlapping by removing it
                            //std::cout << "single subchart overlaps " << id << std::endl;
                            assert(inStack.size() == 0 && s.size() == 0); // we detect this from the starting face
                            for (auto fptr : chfp->fpVec) {
                                visited.insert(fptr); // set all the subchart faces as visited
                            }
                            currentMerge.insert(chfp);
                            break;
                        } else {
                            // if it does, backtrack and forbid the merge
                            aggregate.fpVec.erase(aggregate.fpVec.end() - chfp->fpVec.size(), aggregate.fpVec.end());
                            forbidden.insert(id);
                        }
                    } else {
                        expandVisit = true;
                    }
                }
                if (expandVisit) {
                    visited.insert(fp);
                    if (currentMerge.count(chfp) == 0) {
                        currentMerge.insert(chfp);
                    }

                    for (int i = 0; i < 3; ++i) {
                        if ((inStack.count(fp->FFp(i)) == 0) && (visited.count(fp->FFp(i)) == 0) && (forbidden.count(ccid[fp->FFp(i)]) == 0)) {
                            s.push(fp->FFp(i));
                            inStack.insert(fp->FFp(i));
                        }
                    }
                }
            } // while (!s.empty())
            assert(inStack.empty());
            for (auto cc : currentMerge) {
                for (auto ff : cc->fpVec)
                    assert(visited.count(ff) > 0);
            }
            mergeSets.push_back(currentMerge);
        }
    }

    tri::UpdateTopology<Mesh>::FaceFace(m);

    int k = 0;
    for (auto& set : mergeSets) {
        /*
        std::cout << "Merge set " << k++ << " = ";
        for (auto ch : set) {
            std::cout << ch->id << " ";
        }
        std::cout << endl;
        */
        auto res = gm.Collapse(set.begin(), set.end());
        assert(res.first == GraphManager::Collapse_OK);
        chartVec.push_back(res.second);
    }

    if (binarySplit) {
        while (chartVec.size() > 2) {
            std::sort(chartVec.begin(), chartVec.end(),
                          [](const ChartHandle& c1, const ChartHandle& c2) { return c1->Area3D() < c2->Area3D(); }
            );

            ChartHandle c0 = chartVec[0];
            for (std::size_t i = 1; i < chartVec.size(); ++i) {
                ChartHandle ci = chartVec[i];
                GraphManager::Edge e{c0, chartVec[i]};
                if (gm.ExistsEdge(e)) {
                    ChartHandle merged = gm.Collapse(e);
                    chartVec.erase(std::remove_if(chartVec.begin(), chartVec.end(), [&](const ChartHandle& ch) { return ch == c0 || ch == ci; }),
                            chartVec.end());
                    chartVec.push_back(merged);
                    break;
                }
            }
        }
    }

    /*
    static int n = 0;

    for (std::size_t i = 0; i < chartVec.size(); ++i) {
        std::stringstream ss;
        ss.clear();
        ss << "split_" << n << "_" << i << ".obj";
        Mesh shell;
        BuildShell(shell, *chartVec[i], ParameterizationGeometry::Texture, true);
        MarkInitialSeamsAsFaux(shell, m);
        for (auto& sf : shell.face) {
            for (int i = 0; i < 3; ++i) {
                if (sf.IsF(i)) {
                    sf.C() = vcg::Color4b::Blue;
                    sf.FFp(i)->C() = vcg::Color4b::Blue;
                }
            }
        }
        tri::io::Exporter<Mesh>::Save(shell, ss.str().c_str(), tri::io::Mask::IOM_ALL);
    }
    n++;
    */


}

void PreprocessMesh(Mesh& m, GraphHandle graph)
{
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

    ParameterizationStrategy strategy = DefaultStrategy();
    strategy.directParameterizer = FixedBorderBijective;
    strategy.energy = EnergyType::SymmetricDirichlet;
    strategy.geometry = Model;
    strategy.descent = ScalableLocallyInjectiveMappings;
    strategy.optimizerIterations = 200;
    strategy.padBoundaries = true;
    strategy.applyCut = false;
    strategy.warmStart = false;

    TextureObjectHandle textureObject = graph->textureObject;
    double uvRatio = textureObject->TextureWidth(0) / (double) textureObject->TextureHeight(0);

    for (auto& entry : graph->charts) {
        auto chart = entry.second;
        if (chart->AreaUV() > 0)
            continue;
        numNoParam++;

        std::cout << "Parameterizing region of " << chart->FN() << " zero UV area faces" << std::endl;

        ParameterizerObject po{chart, strategy};
        po.SetEnergyDiffTolerance(0);
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

                    // Handle squashed texture coordinates due to non square textures
                    if (uvRatio > 1)
                        fptr->WT(i).P().X() /= uvRatio;
                    else if (uvRatio < 1)
                        fptr->WT(i).P().Y() *= uvRatio;

                    fptr->WT(i).P().X() -= (box.min.X() + box.DimX());
                    fptr->WT(i).P().Y() -= (box.min.Y() + randomDisplacementV);
                    fptr->WT(i).P() *= scale;
                }
            }
            numParameterized++;
            chart->ParameterizationChanged();
        }
    }

    std::cout << "[LOG] Newly parameterized regions: " << numParameterized << "/" << numNoParam << std::endl;
}

int ParameterizeGraph(GraphManager& gm, ParameterizationStrategy strategy, double injectivityTolerance)
{
    struct ParamTask {
        ChartHandle chart;
        bool warmStart;
    };

    bool injectivityCheckRequired;
    if (strategy.scaffold || injectivityTolerance < 0)
        injectivityCheckRequired = false;
    else {
        assert(injectivityTolerance >= 0 && injectivityTolerance <= 1);
        injectivityCheckRequired = true;
    }

    Timer timer;

    GraphHandle graph = gm.Graph();

    assert(HasWedgeTexCoordStorageAttribute(graph->mesh));
    auto wtcsattr = GetWedgeTexCoordStorageAttribute(graph->mesh);

    double uvRatio = graph->textureObject->TextureWidth(0) / (double) graph->textureObject->TextureHeight(0);

    std::cout << "Parameterizing " << graph->charts.size() << " regions..." << std::endl;

    // gather the charts to parameterize
    std::deque<ParamTask> paramQueue;
    for (auto entry : graph->charts) {
        ChartHandle chart = entry.second;
        if (uvRatio != 1.0 || chart->numMerges > 0)
            // if the texture is non square, put the chart in the queue anyway so that
            // the coordinates are restored with the corrected aspect ratio for the
            // subsequent packing
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
            // restore the initial texture coordinates
            for (auto fptr : chart->fpVec) {
                TexCoordStorage tcs = wtcsattr[fptr];
                for (int i = 0; i < 3; ++i) {
                    fptr->WT(i) = tcs.tc[i];
                }
            }
            continue;
        }

        ParameterizationStrategy strat = strategy;
        if (task.warmStart)
            strat.warmStart = true;

        ParameterizerObject po{chart, strat};
        bool parameterized = po.Parameterize();
        if (parameterized) {
            std::cout << "  Chart parameterized correctly" << std::endl;
            po.SyncChart();

            bool backtrack = false;
            if (injectivityCheckRequired) {
                RasterizedParameterizationStats stats = GetRasterizationStats(chart, 1024, 1024);
                double fraction = stats.lostFragments / (double) stats.totalFragments;
                if (fraction > injectivityTolerance) {
                    std::cout << "WARNING: REGION " << chart->id << " HAS OVERLAPS IN THE PARAMETERIZATION (overlap fraction = " << fraction << ")" << std::endl;
                    backtrack = true;
                }

            }
            if (backtrack) {
                // split the chart in the mesh graph using the graph manager, and readd
                // the new task with the warmStart flag set to true
                std::vector<ChartHandle> splitCharts;
                gm.Split(chart->id, splitCharts);
                std::vector<ChartHandle> newCharts;
                RecoverFromSplit(splitCharts, gm, newCharts, true);
                for (auto& c : newCharts) {
                    paramQueue.push_back(ParamTask{c, strategy.warmStart});
                }
            } else {
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
            }
        } else {
            numFailed++;
            std::cout << "Parameterization failed" << std::endl;
            // split the aggregate and restore the original uvs
            std::vector<ChartHandle> splitCharts;
            gm.Split(chart->id, splitCharts);
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

        std::cout << "Iteration " << iter++ << " took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    }

    std::cout << "Parameterization took " << timer.TimeElapsed() << " seconds" << std::endl;

    return numFailed;
}

void Pack(GraphHandle graph)
{
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
        ChartOutlinesUV(graph->mesh, *chart, uvOutlines);
        int i = tri::OutlineUtil<float>::LargestOutline2(uvOutlines);
        if (tri::OutlineUtil<float>::Outline2Area(uvOutlines[i]) < 0)
            tri::OutlineUtil<float>::ReverseOutline2(uvOutlines[i]);

        outlineMap[chart->id] = texOutlines.size();
        texOutlines.push_back(uvOutlines[i]);
    }


    std::cout << "Packing the atlas..." << std::endl;

    Timer timer;

    // pack the atlas TODO function parameter to choose the packing strategy
    RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters packingParam;
    packingParam.costFunction  = RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters::LowestHorizon;
    //packingParam.costFunction  = RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters::MinWastedSpace;
    //packingParam.costFunction  = RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters::MixedCost;
    packingParam.doubleHorizon = true;
    packingParam.cellSize = 2;
    packingParam.pad = 4;
    packingParam.rotationNum = 16; //number of rasterizations in 90°

    TextureObjectHandle to = graph->textureObject;
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
