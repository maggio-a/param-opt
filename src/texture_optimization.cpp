#include "texture_optimization.h"

#include "uv.h"
#include "mesh_graph.h"
#include "timer.h"
#include "dcpsolver.h"
#include "fixed_border_bijective.h"
#include "iterative.h"
#include "metric.h"
#include "parameterization_checker.h"

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <vcg/space/rasterized_outline2_packer.h>
#include <vcg/space/outline2_packer.h>
#include <wrap/qt/outline2_rasterizer.h>
#include <vcg/space/index/grid_util2d.h>
#include <vcg/space/segment2.h>
#include <vcg/space/intersection2.h>

//#include<vcg/complex/algorithms/cut_tree.h>
//#include<vcg/complex/algorithms/curve_on_manifold.h>
//#include<vcg/complex/algorithms/crease_cut.h>

#include <vector>


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
    ChartOutlinesUV<double>(m, chart, outlines);

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
    // Compute preliminary parameterization graph
    auto graph = ComputeParameterizationGraph(m, nullptr, nullptr);

    // Parameterize regions that are not parameterized
    ReparameterizeZeroAreaRegions(m, graph);

}

void ReparameterizeZeroAreaRegions(Mesh &m, std::shared_ptr<MeshGraph> graph)
{
    double scale;
    DistortionMetric::ComputeAreaScale(m, scale, ParameterizationGeometry::Model);
    scale = std::sqrt(1.0 / scale);

    int numNoParam = 0;
    int numParameterized = 0;

    ParameterizationStrategy strategy;
    strategy.directParameterizer = FixedBorderBijective;
    strategy.optimizer = SymmetricDirichletOpt;
    strategy.geometry = Model;

    for (auto& entry : graph->charts) {
        auto chart = entry.second;

        if (chart->AreaUV() > 0) continue;

        numNoParam++;

        strategy.descent = ScalableLocallyInjectiveMappings;
        strategy.optimizerIterations = 200;
        std::cout << "Parameterizing region of " << chart->FN() << " zero UV area faces" << std::endl;
        bool parameterized = ParameterizeChart(m, chart, strategy);

        if (!parameterized) {
            std::cout << "WARNING: preliminary parameterization of chart " << chart->id << " failed" << std::endl;
        } else {
            /* as convenience to detect regions that originally did not have a parameterization, the texcoords
             * of such regions have negative u and a randomly displaced v */
            Box2d box = chart->UVBox();
            double randomDisplacementU = rand() / (double) RAND_MAX;
            for (auto fptr : chart->fpVec) {
                for (int i = 0; i < 3; ++i) {
                    fptr->WT(i).P().X() -= (box.min.X() + box.DimX());
                    fptr->WT(i).P().Y() -= (box.min.Y() + randomDisplacementU);
                    fptr->WT(i).P() *= scale;
                }
            }
            numParameterized++;
        }
    }

    std::cout << "[LOG] Newly parameterized regions: " << numParameterized << "/" << numNoParam << std::endl;
}

bool ParameterizeMesh(Mesh& m, ParameterizationStrategy strategy)
{
    bool solved = false;
    if (strategy.directParameterizer == DirectParameterizer::DCP) {
        DCPSolver<Mesh> solver(m);
        if (strategy.geometry == ParameterizationGeometry::Model) {
            solved = solver.Solve(DefaultVertexPosition<Mesh>{});
        } else {
            WedgeTexCoordAttributePosition<Mesh> texcoordPosition{m, "WedgeTexCoordStorage"};
            solved = solver.Solve(texcoordPosition);
        }
    } else if (strategy.directParameterizer == DirectParameterizer::FixedBorderBijective) {
        UniformSolver<Mesh> solver(m);
        solved = solver.Solve();
        if (strategy.padBoundaries == false) {
            // delete the faces that were added by the hole filling process
            for (auto & f : m.face) {
                if (f.IsS()) tri::Allocator<Mesh>::DeleteFace(m, f);
            }
            tri::Allocator<Mesh>::CompactEveryVector(m);
        }
        // check for inverted faces
        //for (auto& f : m.face) {
        //    double areaUV = (f.V(1)->T().P() - f.V(0)->T().P()) ^ (f.V(2)->T().P() - f.V(0)->T().P());
        //    assert(areaUV > 0 && "Parameterization is not bijective");
        //}
    } else {
        assert(0 && "Solver not supported");
    }



    if (solved) { // optimize and Copy texture coords back

        if (strategy.optimizerIterations > 0) {

            //GradientDescent opt(std::make_shared<SymmetricDirichlet>(pm));
            Energy::Geometry mode = Energy::Geometry::Model;
            if (strategy.geometry == ParameterizationGeometry::Texture) mode = Energy::Geometry::Texture;

            DescentMethod *opt;
            switch(strategy.descent) {
            case DescentType::Gradient:
                opt = new GradientDescent{std::make_shared<SymmetricDirichlet>(m, mode)};
                break;
            case DescentType::LimitedMemoryBFGS:
                opt = new LBFGS{std::make_shared<SymmetricDirichlet>(m, mode), 10};
                break;
            case DescentType::ScalableLocallyInjectiveMappings:
                opt = new SLIM{std::make_shared<SymmetricDirichlet>(m, mode)};
                break;
            default:
                assert(0);
            }

            Timer t;
            int i;
            double energyVal, gradientNorm, energyDiff;
            for (i = 0; i < strategy.optimizerIterations; ++i) {
                energyVal = opt->Iterate(gradientNorm, energyDiff);
                for (auto& f : m.face) {
                    double areaUV = (f.V(1)->T().P() - f.V(0)->T().P()) ^ (f.V(2)->T().P() - f.V(0)->T().P());
                    assert(areaUV > 0 && "Parameterization is not bijective");
                }
                if (gradientNorm < 1e-3) {
                    std::cout << "Stopping because gradient magnitude is small enough (" << gradientNorm << ")" << std::endl;
                    break;
                }
                if (energyDiff < 1e-9) {
                    std::cout << "Stopping because energy improvement is too small (" << energyDiff << ")" << std::endl;
                    break;
                }
            }
            if (i >= strategy.optimizerIterations) {
                std::cout << "Maximum number of iterations reached" << std::endl;
            }
            std::cout << "Stopped after " << i << " iterations, gradient magnitude = " << gradientNorm
                      << ", energy value = " << energyVal << std::endl;

            std::cout << "Optimization took " << t.TimeSinceLastCheck() << " seconds" << std::endl;

            tri::UpdateTexture<Mesh>::WedgeTexFromVertexTex(m);

            delete opt;

        }
    }
    return solved;
}

bool ParameterizeChart(Mesh &m, GraphManager::ChartHandle ch, ParameterizationStrategy strategy)
{
    Mesh pm;
    bool sanitize = (strategy.directParameterizer == DirectParameterizer::FixedBorderBijective);
    bool copyWedgeTexCoordStorage = (strategy.geometry == ParameterizationGeometry::Texture);
    CopyFaceGroupIntoMesh(pm, *ch, sanitize, copyWedgeTexCoordStorage);

    bool solved = ParameterizeMesh(pm, strategy);
    if (solved) {

        for (std::size_t i = 0; i < ch->fpVec.size(); ++i) {
            for (int k = 0; k < 3; ++k) {
                ch->fpVec[i]->WT(k).P() = pm.face[i].WT(k).P();
            }
        }
    }
    return solved;
}

/* Ugly function that tries to perform subsequent merges after a split. split is the vector of charts
 * that formed the aggregate, chartQueue is the queue were the newly merged charts must be inserted
 * */
static void RecoverFromSplit(std::vector<ChartHandle>& split, GraphManager& gm, std::deque<ChartHandle>& chartQueue)
{
    /* very simple heuristic, select the two largest charts in the split, and iteratively grow them
     * until all the charts in the split are covered, producing to 2 new aggregates that will be
     * parameterized independently */
    std::sort(split.begin(), split.end(),
              [](const ChartHandle& c1, const ChartHandle& c2) { return c1->Area3D() > c2->Area3D(); }
    );
    ChartHandle c1 = split[0];
    ChartHandle c2 = split[1];
    std::unordered_set<ChartHandle> charts;
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

    chartQueue.push_back(c1);
    chartQueue.push_back(c2);
}

int ParameterizeGraph(GraphManager& gm,
                      ParameterizationStrategy strategy,
                      bool failsafe,
                      double threshold,
                      bool retry)
{
    assert(threshold >= 0 && threshold <= 1);

    Timer timer;

    auto graph = gm.Graph();
    Mesh& m = graph->mesh;

    auto WTCSh = tri::Allocator<Mesh>::FindPerFaceAttribute<TexCoordStorage>(m, "WedgeTexCoordStorage");
    assert(tri::Allocator<Mesh>::IsValidHandle<TexCoordStorage>(m, WTCSh));

    std::cout << "Parameterizing " << graph->charts.size() << " regions..." << std::endl;

    // gather the charts to parameterize
    std::deque<ChartHandle> chartQueue;
    for (auto entry : graph->charts) {
        ChartHandle chart = entry.second;
        if (chart->numMerges > 0) chartQueue.push_back(chart);
        chart->ParameterizationChanged(); // packing changes texture coords even for charts that are not reparameterized

        /* !!!
         * If the chart was not merged because it was a disconnected component AND a closed surface, the outline
         * extraction will fail (this means that the original parameterization mapped a closed surface to a flat area)
         * a possible workaround would be to use the uv bounding box instead of the outline for those cases
         * */
    }

    int iter = 0;
    int numFailed = 0;
    while (chartQueue.size() > 0) {
        ChartHandle chart = chartQueue.front();
        chartQueue.pop_front();

        std::cout << "Chart " << chart->id << " - FN=" << chart->FN() << ", FI=" << tri::Index(m, chart->Fp()) << std::endl;

        if (chart->numMerges == 0) continue;

        double oldUvArea = chart->AreaUV();

        bool needToSplit = false;
        bool parameterized = ParameterizeChart(m, chart, strategy);
        if (!parameterized) {
            numFailed++;
            std::cout << "Parameterization failed" << std::endl;
            needToSplit = true; // split if parameterization failed
        }
        else {
            if (failsafe) {
                RasterizedParameterizationStats stats = GetRasterizationStats(chart, 1024, 1024);
                double fraction = stats.lostFragments / (double) stats.totalFragments;
                if (fraction > threshold) {
                    std::cout << "WARNING: REGION " << chart->id << " HAS OVERLAPS IN THE PARAMETERIZATION "
                              << "(overlap fraction = " << fraction << ")" << std::endl;
                    needToSplit = true; // split if overlaps are above threshold
                }
            }
        }

        if (needToSplit) {
            // restore original texture coordinates for each face of the chart
            for (auto fptr : chart->fpVec) {
                TexCoordStorage tcs = WTCSh[fptr];
                for (int i = 0; i < 3; ++i) {
                    fptr->WT(i) = tcs.tc[i];
                    fptr->V(i)->T() = tcs.tc[i];
                }
            }
            // split the chart in the mesh graph using the graph manager
            std::vector<ChartHandle> splitCharts;
            gm.Split(chart->id, splitCharts);
            if (retry) {
                RecoverFromSplit(splitCharts, gm, chartQueue);
            }
        }
        else {
            // if the parameterization is valid, scale chart relative to the original parameterization area value

            // Normalize area: the region gets scaled by sqrt(oldUVArea)/sqrt(newUVArea) to keep the original proportions
            // between the regions of the atlas in an attempt to not lose too much texture data when the new texture is rendered
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

        std::cout << "Iteration " << iter++ << " took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    }

    // Pack the atlas
    std::vector<std::vector<Point2f>> texOutlines;
    std::unordered_map<RegionID,std::size_t> outlineMap; // map each region to the index of its outline in texOutlines

    for (auto entry : graph->charts) {
        GraphManager::ChartHandle chart = entry.second;
        std::cout << "Chart " << chart->id << " - FN=" << chart->FN() << ", FI=" << tri::Index(m, chart->Fp()) << std::endl;

        // Save the outline of the parameterization for this portion of the mesh
        std::vector<std::vector<Point2f>> uvOutlines;
        ChartOutlinesUV(m, chart, uvOutlines);
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
    packingParam.cellSize = 4;
    packingParam.rotationNum = 16; //number of rasterizations in 90°

    Point2i gridSize(1024, 1024);
    std::vector<Similarity2f> transforms;

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
                fptr->WT(j).P() = Point2d{transformedTexCoordPos[0] / 1024.0, transformedTexCoordPos[1] / 1024};
            }
        }
    }

    return numFailed;
}


void ReduceTextureFragmentation(Mesh &m, std::shared_ptr<MeshGraph> graph, std::size_t minRegionSize)
{
    if (minRegionSize == 0) return;

    assert(minRegionSize < (std::size_t) m.FN());

    tri::UpdateTopology<Mesh>::FaceFace(m);

    GraphManager gm{graph};

    std::cout << "Initialized graph manager" << std::endl;

    ReduceTextureFragmentation_NoPacking(gm, minRegionSize);

    assert(0 && "TODO ParameterizeGraph parameters");
    int c = ParameterizeGraph(gm, ParameterizationStrategy{}, true, 0);
    if (c > 0) std::cout << "WARNING: " << c << " regions were not parameterized correctly" << std::endl;
}


void ReduceTextureFragmentation_NoPacking(GraphManager &gm, std::size_t minRegionSize)
{
    Timer timer;
    int mergeCount;
    int numIter = 0;

    do {
        mergeCount = gm.CloseMacroRegions(minRegionSize);

        while (gm.HasNextEdge()) {
            auto we = gm.PeekNextEdge();
            if (we.first.a->FN() > minRegionSize && we.first.b->FN() > minRegionSize)
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
