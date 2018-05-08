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
#include <vcg/space/index/grid_util2d.h>
#include <vcg/space/segment2.h>
#include <vcg/space/intersection2.h>

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
    std::cout << "FIXME" << std::endl;
    // Compute preliminary parameterization graph
    auto graph = ComputeParameterizationGraph(m, nullptr);

    // Parameterize regions that are not parameterized
    ReparameterizeZeroAreaRegions(m, graph);

}

void ReparameterizeZeroAreaRegions(Mesh &m, std::shared_ptr<MeshGraph> graph)
{

    std::cout << "TODO FIXME" << std::endl;
    return;

#if 0
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
        bool parameterized = ParameterizeChart(m, chart, strategy, false);

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
#endif
}

struct EdgeDistortionCounter {
    Mesh::FacePointer fp;
    int i;
    double distortion;

    EdgeDistortionCounter(Mesh::FacePointer fptr, int ii) : fp{fptr}, i{ii}, distortion{0}
    {
    }
};


struct FacePair {
    int i1;
    int i2;

    bool operator==(const FacePair& other) const
    {
        return (i1 == other.i1 && i2 == other.i2) || (i1 == other.i2 && i2 == other.i1);
    }
};

struct IntPairHasher {
    std::size_t operator()(const std::pair<int, int>& p) const noexcept
    {
        std::size_t seed = 0;
        int a = std::min(p.first, p.second);
        int b = std::max(p.first, p.second);
        seed ^= std::hash<int>()(a) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<int>()(b) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

#include <wrap/io_trimesh/export.h>

bool ParameterizeShell(Mesh& shell, ParameterizationStrategy strategy, Mesh& baseMesh)
{
    if (strategy.warmStart) {
        // if warm start is requested, make sure that the initial parameterization is injective
        // (should also check that it is connected...)
        assert(CheckLocalInjectivity(shell));
        assert(CheckUVConnectivity(shell));
    } else {
        // compute initial bijective parameterization (Tutte)
        UniformSolver<Mesh> solver(shell);
        bool solved = solver.Solve();
        if (!solved)
            return false;
    }

    if (strategy.padBoundaries == false)
        ClearHoleFillingFaces(shell);

    if (strategy.optimizerIterations > 0) {
        Timer t;
        int i;
        double energyVal, normalizedEnergyVal, gradientNorm, energyDiff;
        std::cout << "TODO optimize loop ParameterizeShell(): avoid recreating descent and energy objects" << std::endl;

        bool appliedCut = false;
        for (i = 0; i < strategy.optimizerIterations; ++i) {
            tri::io::Exporter<Mesh>::Save(shell, "shell.obj", tri::io::Mask::IOM_FACECOLOR | tri::io::Mask::IOM_VERTTEXCOORD);

            // create energy and optimizer
            auto energy = std::make_shared<SymmetricDirichlet>(shell);
            std::shared_ptr<DescentMethod> opt;
            switch(strategy.descent) {
            case DescentType::Gradient:
                opt = std::make_shared<GradientDescent>(energy);
                break;
            case DescentType::LimitedMemoryBFGS:
                opt = std::make_shared<LBFGS>(energy, 10);
                break;
            case DescentType::ScalableLocallyInjectiveMappings:
                opt = std::make_shared<SLIM>(energy);
                break;
            default:
                assert(0);
            }

            energyVal = opt->Iterate(gradientNorm, energyDiff);

            for (auto& sf : shell.face) {
                double areaUV = (sf.V(1)->T().P() - sf.V(0)->T().P()) ^ (sf.V(2)->T().P() - sf.V(0)->T().P());
                assert(areaUV > 0 && "Parameterization is not bijective");
            }
            normalizedEnergyVal = energy->E_IgnoreMarkedFaces(true);
            if (gradientNorm < 1e-3) {
                std::cout << "Stopping because gradient magnitude is small enough (" << gradientNorm << ")" << std::endl;
                break;
            }
            if (energyDiff < 1e-9) {
                std::cout << "Stopping because energy improvement is too small (" << energyDiff << ")" << std::endl;
                break;
            }

            tri::UpdateColor<Mesh>::PerFaceQualityRamp(shell);
            SyncShell(shell);
            if (i > 0 && (i % 10) == 0) {
                RemeshShellHoles(shell, strategy.geometry, baseMesh);
                energy->UpdateCache();
            }

            tri::UpdateTexture<Mesh>::WedgeTexFromVertexTex(shell);
            energy->MapToFaceQuality(true);


            if (i > 25 && i > (strategy.optimizerIterations / 2) && appliedCut == false) {
                tri::UpdateTexture<Mesh>::WedgeTexFromVertexTex(shell);
                energy->UpdateCache();
                energy->MapToFaceQuality(true);

                // Mark texture seams in the shell mesh along which the procedure is
                // allowed to cut in order to relieve the distortion of the flattening
                // if it goes above a given threshold. Note that if a face is
                // hole-Filling, the energy computation automatically attenuates distortion
                MarkInitialSeamsAsFaux(shell, baseMesh);

                // code to mark all edges as candidates for the shortest path
                //tri::UpdateTopology<Mesh>::FaceFace(shell);
                //tri::UpdateFlags<Mesh>::FaceBorderFromFF(shell);
                //tri::UpdateFlags<Mesh>::FaceSetF(shell);
                //for (auto& sf : shell.face) {
                //    for (int i = 0; i < 3; ++i) {
                //        if (sf.IsB(i)) sf.ClearF(i);
                //    }
                //}
                double maxDistance = ComputeDistanceFromBorderOnSeams(shell);

                tri::UpdateColor<Mesh>::PerFaceQualityRamp(shell);
                float minq, maxq;
                tri::Stat<Mesh>::ComputePerFaceQualityMinMax(shell, minq, maxq);
                std::cout << "Min distortion value = " << minq << std::endl;
                std::cout << "Max distortion value = " << maxq << std::endl;

                // Select candidate crease edge for cutting

                double maxEnergy = 0;
                Mesh::FacePointer startFace = nullptr;
                int cutEdge = -1;
                for (auto& sf : shell.face) {
                    for (int i = 0; i < 3; ++i) if (sf.IsF(i)) {
                        double w = std::max(sf.V0(i)->Q(), sf.V1(i)->Q()) / maxDistance;
                        double weightedEnergy = energy->E(sf, true) * w;
                        // w is INFINITY if the seam does not reach the mesh boundary
                        if (std::isfinite(weightedEnergy) && weightedEnergy > maxEnergy) {
                            maxEnergy = weightedEnergy;
                            startFace = &sf;
                            cutEdge = i;
                        }
                    }
                }
                assert(startFace != nullptr);

                PosF p(startFace, cutEdge);
                assert(!p.IsBorder());

                SelectShortestSeamPathToBoundary(shell, p);

                tri::CutMeshAlongSelectedFaceEdges(shell);
                tri::UpdateTopology<Mesh>::FaceFace(shell);
                tri::UpdateFlags<Mesh>::VertexBorderFromFaceAdj(shell);
                for (auto& sf : shell.face) {
                    for (int i = 0; i < 3; i++) {
                        if (face::IsBorder(sf, i)) sf.ClearF(i);
                    }
                }

                tri::io::Exporter<Mesh>::Save(shell, "shell_cut.obj", tri::io::Mask::IOM_FACECOLOR);
                appliedCut = true;
            }
        }
        std::cout << "Stopped after " << i << " iterations, gradient magnitude = " << gradientNorm
                  << ", normalized energy value = " << normalizedEnergyVal << std::endl;

        float minq, maxq;
        tri::Stat<Mesh>::ComputePerFaceQualityMinMax(shell, minq, maxq);
        std::cout << "Min distortion value = " << minq << std::endl;
        std::cout << "Max distortion value = " << maxq << std::endl;

        std::cout << "Optimization took " << t.TimeSinceLastCheck() << " seconds" << std::endl;
    }

    return true;
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
    std::cout << "Recovery produced two charts of sizes " << c1->numMerges + 1 << " " << c2->numMerges + 1 << std::endl;
}

/*
 * if overlap detected
 * split chart
 * recover from split
 * for each filled hole
 *   if the hole crosses the split, remove it
 * for each edge
 *   if the edge lies on the split, duplicate it
 * the pm is now made of two distinct connected components, detach them and parameterize each independently
 * using the previously computed solution as starting point
 */
bool ParameterizeChart(Mesh &m, ChartHandle ch, ParameterizationStrategy strategy)
{
    Mesh shell;
    return ParameterizeChart(m, ch, strategy, shell);
}

bool ParameterizeChart(Mesh& m, ChartHandle ch, ParameterizationStrategy strategy, Mesh& shell)
{
    shell.Clear();
    BuildShell(shell, *ch, strategy.geometry);

    std::cout << "WARNING forcing warm start to true" << std::endl;
    strategy.warmStart = true;

    bool solved = ParameterizeShell(shell, strategy, ch->mesh);
    if (solved) {
        for (std::size_t i = 0; i < ch->fpVec.size(); ++i) {
            for (int k = 0; k < 3; ++k) {
                ch->fpVec[i]->WT(k).P() = shell.face[i].WT(k).P();
            }
        }
    }
    return solved;

}

int ParameterizeGraph(GraphManager& gm,
                      double packingCoverage,
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

        std::cout << "+ Chart parameterized correctly" << std::endl;

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
            if (strategy.geometry == Model || strategy.optimizerIterations == 0) {
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

        chart->ParameterizationChanged(); // packing changes texture coords even for charts that are not reparameterized

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
    packingParam.cellSize = 2;
    packingParam.pad = 16;
    //packingParam.pad = 0;
    packingParam.rotationNum = 16; //number of rasterizations in 90°

    TextureObjectHandle to = gm.Graph()->textureObject;
    double scale = std::sqrt(packingCoverage);
    int gridWidth = to->TextureWidth(0) / 2;
    int gridHeight = to->TextureHeight(0) / 2;
    Point2i gridSize(gridWidth, gridHeight);
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
                fptr->WT(j).P() = Point2d{transformedTexCoordPos[0] / double(gridSize.X()), transformedTexCoordPos[1] / double(gridSize.Y())};
                fptr->WT(j).P() *= scale;
            }
        }
    }

    return numFailed;
}


/*
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
*/

void ReduceTextureFragmentation_NoPacking_TargetRegionCount(GraphManager &gm, std::size_t regionCount, std::size_t minRegionSize)
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
            bool regionReached = gm.Graph()->Count() <= regionCount;
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


void ReduceTextureFragmentation_NoPacking(GraphManager &gm, std::size_t minRegionSize)
{
    Timer timer;
    int mergeCount;
    int numIter = 0;

    std::cout << "[LOG] Reduction strategy MinRegionSize=" << minRegionSize << std::endl;

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
