#include "texture_optimization.h"

#include "mesh.h"
#include "uv.h"
#include "mesh_graph.h"
#include "timer.h"
#include "iterative_solvers.h"
#include "metric.h"
#include "mesh_utils.h"
#include "packing_utils.h"
#include "texture_rendering.h"
#include "logging.h"

#include "parallel.h"

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

#include <vcg/complex/algorithms/hole.h>
#include <vcg/complex/algorithms/isotropic_remeshing.h>

#include <wrap/io_trimesh/export.h>

#include <vector>
#include <algorithm>

#include <vcg/complex/algorithms/update/color.h>


void RemoveDegeneracies(Mesh& m)
{
    int dupVert = tri::Clean<Mesh>::RemoveDuplicateVertex(m);
    //int zeroArea = tri::Clean<Mesh>::RemoveZeroAreaFace(m);

    tri::Allocator<Mesh>::CompactEveryVector(m);

    tri::UpdateTopology<Mesh>::FaceFace(m);

    if (dupVert > 0)
        LOG_VERBOSE << "Removed " << dupVert << " duplicate vertices";

    //if (zeroArea > 0)
    //    LOG_VERBOSE << "Removed " << zeroArea << " zero area faces";

    int numVertexSplit = 0;
    int nv;
    while ((nv = tri::Clean<Mesh>::SplitNonManifoldVertex(m, 0)) > 0)
        numVertexSplit += nv;
    if (numVertexSplit > 0)
        LOG_VERBOSE << "Mesh was not vertex manifold, split " << numVertexSplit << " vertices";
    else
        LOG_VERBOSE << "Mesh is vertex manifold";

    int numRemovedFaces = tri::Clean<Mesh>::RemoveNonManifoldFace(m);
    if (numRemovedFaces > 0)
        LOG_VERBOSE << "Mesh was not edge manifold, removed " << numRemovedFaces << " faces";
    else
        LOG_VERBOSE << "Mesh is edge manifold";

    tri::Allocator<Mesh>::CompactEveryVector(m);

    // Handle zero area faces by selecting them, dilating the selation and remeshing
    // the selected areas. Zero area holes are handles in the same way by being filled
    // beforehand

    LOG_INFO << "FIXME Remeshing does not guarantee that texture coordinates are preserved";
    tri::UpdateSelection<Mesh>::Clear(m);

    unsigned fn = m.face.size();
    tri::Hole<Mesh>::EarCuttingFill<tri::MinimumWeightEar<Mesh>>(m, 4, false);

    // If there are isolated tris, the hole-filling algorithm adds degenerate faces, so remove them
    // before iterating
    tri::Clean<Mesh>::RemoveDegenerateFace(m);

    for (unsigned i = 0; i < m.face.size(); ++i) {
        if (!m.face[i].IsD()) {
            double farea = DoubleArea(m.face[i]) / 2.0;
            if ((farea == 0) || ((i >= fn) && (farea < 1e-8))) {
                m.face[i].SetS();
            } else if (i >= fn) {
                LOG_DEBUG << "Deleting hole face of area " << farea;
                tri::Allocator<Mesh>::DeleteFace(m, m.face[i]);
            }
        }
    }
    tri::Clean<Mesh>::RemoveUnreferencedVertex(m);
    tri::Allocator<Mesh>::CompactEveryVector(m);

    // dilate selection
    tri::UpdateSelection<Mesh>::VertexFromFaceLoose(m, true);
    tri::UpdateSelection<Mesh>::FaceFromVertexLoose(m, true);

    double avgEdgeLen = 0;
    int selCount = 0;
    for (auto& f : m.face) {
        if (f.IsS()) {
            for (int i = 0; i < 3; ++i) {
                avgEdgeLen += EdgeLength(f, i);
                selCount++;
            }
        }
    }
    avgEdgeLen /= selCount;

    if (selCount > 0) {
        LOG_DEBUG << "Detected zero area faces/holes";
        IsotropicRemeshing<Mesh>::Params params;
        params.stat.Reset();
        params.SetTargetLen(avgEdgeLen);
        params.SetFeatureAngleDeg(30);
        params.splitFlag = true;
        params.swapFlag = true;
        params.collapseFlag = true;
        params.selectedOnly = true;
        params.smoothFlag = false;
        params.projectFlag = false;
        params.iter = 1;
        IsotropicRemeshing<Mesh>::Do(m, params);
        tri::Allocator<Mesh>::CompactEveryVector(m);
        tri::UpdateTopology<Mesh>::FaceFace(m);
        LOG_DEBUG << "Zero area face remeshing: " << params.stat.collapseNum << "c " << params.stat.flipNum << "f " << params.stat.splitNum << "s";
    }

    LOG_DEBUG << "Cleaning done";
}

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
                    LOG_DEBUG << "( " <<  a1[0] << " , " << a1[1] << " ) (" << b1[0] << " , " << b1[1] << " )";
                    LOG_DEBUG << "( " <<  a2[0] << " , " << a2[1] << " ) (" << b2[0] << " , " << b2[1] << " )";
                    LOG_DEBUG << intersectionPoint[0] << " , " << intersectionPoint[1];
                    return true;
                }
            }
        }
    }

    return false;
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
            while (!s.empty()) {
                Mesh::FacePointer fp = s.top();
                s.pop();
                ensure_condition(visited.count(fp) == 0);

                ensure_condition(inStack.count(fp) > 0);
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
                            ensure_condition(inStack.size() == 0 && s.size() == 0); // we detect this from the starting face
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
            ensure_condition(inStack.empty());
            for (auto cc : currentMerge) {
                for (auto ff : cc->fpVec)
                    ensure_condition(visited.count(ff) > 0);
            }
            mergeSets.push_back(currentMerge);
        }
    }

    tri::UpdateTopology<Mesh>::FaceFace(m);

    for (auto& set : mergeSets) {
        auto res = gm.Collapse(set.begin(), set.end());
        ensure_condition(res.first == GraphManager::Collapse_OK);
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
}

void RecoverFromFailedInit(std::vector<ChartHandle>& split, GraphManager& gm, std::vector<ChartHandle>& chartQueue)
{
    /* very simple heuristic, select the two largest charts in the split, and iteratively grow them
     * until all the charts in the split are covered, producing to 2 new aggregates that will be
     * parameterized independently */
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

        bool merged = false;

        {
            // grow c1
            ChartHandle bestCollapseCandidate = nullptr;
            double minWeight = std::numeric_limits<double>::infinity();

            for (auto ch : c1->adj) {
                if (charts.count(ch) == 1) {
                    GraphManager::Edge e(c1, ch);
                    double w = gm.EdgeWeight(e);
                    if (w < minWeight) {
                        minWeight = w;
                        bestCollapseCandidate = ch;
                    }
                }
            }
            if (bestCollapseCandidate != nullptr) {
                charts.erase(bestCollapseCandidate);
                c1 = gm.Collapse(GraphManager::Edge(c1, bestCollapseCandidate));
                merged = true;
            }
        }

        {
        // grow c2
            ChartHandle bestCollapseCandidate = nullptr;
            double minWeight = std::numeric_limits<double>::infinity();

            for (auto ch : c2->adj) {
                if (charts.count(ch) == 1) {
                    GraphManager::Edge e(c2, ch);
                    double w = gm.EdgeWeight(e);
                    if (w < minWeight) {
                        minWeight = w;
                        bestCollapseCandidate = ch;
                    }
                }
            }
            if (bestCollapseCandidate != nullptr) {
                charts.erase(bestCollapseCandidate);
                c2 = gm.Collapse(GraphManager::Edge(c2, bestCollapseCandidate));
                merged = true;
            }
        }

        ensure_condition(merged);
    }

    chartQueue.push_back(c1);
    chartQueue.push_back(c2);
}

void ParameterizeZeroAreaRegions(Mesh &m, std::shared_ptr<MeshGraph> graph)
{
    ensure_condition(HasParameterizationScaleInfoAttribute(m));
    auto info = GetParameterizationScaleInfoAttribute(m);
    double scale = info().scale;

    int numNoParam = 0;
    int numParameterized = 0;

    ParameterizationStrategy strategy_cm = MakeStrategy(
            DirectParameterizer::FixedBorderBijective,
            EnergyType::SymmetricDirichlet,
            ParameterizationGeometry::Model,
            DescentType::CompositeMajorization,
            200,            // Number of iterations
            true,           // Fill holes
            false,          // Do not use cuts
            false,          // No warm start
            false           // Disable scaffold
    );

    ParameterizationStrategy strategy_slim = MakeStrategy(
            DirectParameterizer::FixedBorderBijective,
            EnergyType::SymmetricDirichlet,
            ParameterizationGeometry::Model,
            DescentType::ScalableLocallyInjectiveMappings,
            200,            // Number of iterations
            true,           // Fill holes
            false,           // Do not use cuts
            false,          // No warm start
            false           // Disable scaffold
    );


    for (auto& entry : graph->charts) {
        auto chart = entry.second;
        if (chart->AreaUV() > 0)
            continue;
        numNoParam++;

        LOG_DEBUG << "Parameterizing region of " << chart->FN() << " zero UV area faces";

        bool parameterized = false;

        if (chart->FN() == 1) {
            auto fptr = chart->fpVec[0];
            Point2d v10, v20;
            LocalIsometry(fptr->P(1) - fptr->P(0), fptr->P(2) - fptr->P(0), v10, v20);

            fptr->WT(0).P() = Point2d::Zero();
            fptr->WT(1).P() = Point2d::Zero() + v10;
            fptr->WT(2).P() = Point2d::Zero() + v20;

            parameterized = true;
        } else {
            {
                ParameterizerObject po{chart, strategy_cm};
                po.Initialize();
                if (po.GetStatus() == ParameterizerObject::Status::Initialized)
                    parameterized = po.Parameterize();
                if (parameterized)
                    po.SyncChart();
            }

            if (!parameterized) {
                LOG_DEBUG << "Parameterization using CM failed, falling back to SLIM";
                ParameterizerObject po{chart, strategy_slim};
                po.Initialize();
                if (po.GetStatus() == ParameterizerObject::Status::Initialized)
                    parameterized = po.Parameterize();
                if (parameterized)
                    po.SyncChart();
            }
        }

        if (!parameterized) {
            LOG_ERR << "Preliminary parameterization of chart " << chart->id << " failed";
            ensure_condition(0 && "Failed to assign valid parameterization to null chart");
        } else {
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
            chart->ParameterizationChanged();
        }
    }

    LOG_INFO << "Computed new parameterizations for " << numParameterized << " null charts";
}

int RemoveOutliers(GraphHandle& graph)
{
    double totalArea3D = graph->Area3D();
    int count = 0;
    for (auto & entry : graph->charts) {
        ChartHandle c = entry.second;

        double frac = c->Area3D() / totalArea3D;

        if ((frac < 0.0001) && (c->AreaUV() > 0.1)) {
            for (auto fptr : c->fpVec) {
                for (int i = 0; i < 3; ++i) {
                    fptr->WT(i).P() = Point2d(0, 0);
                }
            }
            count++;
        }
    }

    graph = nullptr;

    return count;
}

void RecomputeSegmentation(GraphManager &gm, std::size_t regionCount, double smallRegionAreaThreshold)
{
    Timer timer;
    int mergeCount;
    int numIter = 0;

    ensure_condition(regionCount > 0);

    LOG_INFO << "Clustering: TargetRegionCount=" << regionCount << " , regionThreshold=" << smallRegionAreaThreshold;

    double minChartArea = smallRegionAreaThreshold * gm.Graph()->Area3D();

    do {
        mergeCount = gm.CloseMacroRegions(smallRegionAreaThreshold);

        while (gm.HasNextEdge()) {
            auto we = gm.PeekNextEdge();
            double minNextArea = std::min(we.first.a->Area3D(), we.first.b->Area3D());
            bool regionReached = (regionCount <= 0) || (gm.Graph()->Count() <= regionCount);
            bool sizeThresholdReached = minNextArea > minChartArea;
            if (regionReached && sizeThresholdReached) {
                break;
            } else {
                gm.RemoveNextEdge();
                gm.Collapse(we.first);
                mergeCount++;
                if (mergeCount%50 == 0) {
                    LOG_VERBOSE << "Merged " << mergeCount << " regions...";
                }
            }
        }

        LOG_INFO << "Iteration "  << numIter << " took " << timer.TimeSinceLastCheck()
                 << " seconds (" << mergeCount << " merges took place)";

        numIter++;

    } while (mergeCount > 0);

    // Make sure the texcoords of each chart refer to the same texture unit,
    // otherwise it causes trouble when computing FaceFaceFromTexCoord adjancency
    for (auto entry : gm.Graph()->charts) {
        ChartHandle chart = entry.second;
        for (auto fptr : chart->fpVec)
            for (int i = 0; i < 3; ++i)
                fptr->WT(i).N() = 0;
    }

    LOG_INFO << "Stopping after " << numIter << " passes and " << timer.TimeElapsed() << " seconds";
}

int ParameterizeGraph(GraphManager& gm, ParameterizationStrategy strategy, double injectivityTolerance, int numWorkers)
{
    parallel::Init();

    Timer t;
    parallel::WorkerPool pool(numWorkers, gm);
    pool.Run(strategy, injectivityTolerance);
    LOG_INFO << "Multithreaded parameterization took " << t.TimeElapsed() << " seconds";
    return 0;
}

#if 0
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
        ensure_condition(injectivityTolerance >= 0 && injectivityTolerance <= 1);
        injectivityCheckRequired = true;
    }

    Timer timer;

    GraphHandle graph = gm.Graph();

    ensure_condition(HasWedgeTexCoordStorageAttribute(graph->mesh));
    auto wtcsattr = GetWedgeTexCoordStorageAttribute(graph->mesh);

    std::cout << "Parameterizing " << graph->charts.size() << " regions..." << std::endl;

    // gather the charts to parameterize
    std::deque<ParamTask> paramQueue;
    for (auto entry : graph->charts) {
        ChartHandle chart = entry.second;
        /*
        Mesh probe;
        MeshFromFacePointers(chart->fpVec, probe);
        if (Parameterizable(probe))
            paramQueue.push_back(ParamTask{chart, false});
        */
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

        bool parameterized = false;

        ParameterizerObject po{chart, strat};
        po.Initialize();
        if (po.GetStatus() == ParameterizerObject::Status::Initialized)
            parameterized = po.Parameterize();

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
                ensure_condition(scale > 0);
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
            std::cout << "Parameterization failed" << std::endl;
            // split the aggregate and restore the original uvs
            bool recover = (chart->numMerges > 0);

            std::cout << "Splitting" << std::endl;
            std::vector<ChartHandle> splitCharts;
            gm.Split(chart->id, splitCharts);

            std::cout << "Done" << std::endl;
            if (recover) {
                std::cout << "Recovering" << std::endl;
                std::vector<ChartHandle> newCharts;
                RecoverFromFailedInit(splitCharts, gm, newCharts);
                std::cout << "Done" << std::endl;
                for (auto& c : newCharts) {
                    paramQueue.push_back(ParamTask{c, false});
                }
            } else {
                std::cout << "Chart cannot be split, restoring uvs..." << std::endl;
                numFailed++;
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

        std::cout << "Chart " << chart->id << " parameterization took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    }


    std::cout << "Graph parameterization took " << timer.TimeElapsed() << " seconds" << std::endl;

    return numFailed;
}
#endif

static std::vector<double> ComputeContainerRatios(TextureObjectHandle textureObject, bool oneContainer)
{
    int nTex = oneContainer ? 1 : textureObject->ArraySize();
    std::vector<TextureSize> texSizeVec = ComputeSizes(nTex, textureObject);
    std::vector<double> containerRatios;
    for (auto& tsz : texSizeVec) {
        containerRatios.push_back(tsz.w / (double)tsz.h);
    }
    return containerRatios;
}

RasterizationBasedPacker::PackingStats Pack(GraphHandle graph, const PackingOptions& options)
{
    // Pack the atlas
    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(graph->mesh);

    std::vector<Outline2f> texOutlines;
    std::unordered_map<RegionID,std::size_t> outlineMap; // map each region to the index of its outline in texOutlines

    LOG_VERBOSE << "Computing chart outlines...";

    for (auto entry : graph->charts) {
        GraphManager::ChartHandle chart = entry.second;
        LOG_DEBUG << "Chart " << chart->id << " - FN=" << chart->FN() << ", FI=" << tri::Index(graph->mesh, chart->Fp());

        // Save the outline of the parameterization for this portion of the mesh
        std::vector<Outline2f> uvOutlines;
        ChartOutlinesUV(graph->mesh, *chart, uvOutlines);
        int i = tri::OutlineUtil<float>::LargestOutline2(uvOutlines);
        if (tri::OutlineUtil<float>::Outline2Area(uvOutlines[i]) < 0)
            tri::OutlineUtil<float>::ReverseOutline2(uvOutlines[i]);

        outlineMap[chart->id] = texOutlines.size();
        texOutlines.push_back(uvOutlines[i]);
    }

    //tri::UpdateTopology<Mesh>::FaceFace(graph->mesh);

    LOG_VERBOSE << "Packing charts...";

    Timer timer;

    RasterizationBasedPacker::Parameters packingParam;
    packingParam.costFunction = options.costFunction;
    packingParam.doubleHorizon = true;
    packingParam.innerHorizon = options.useInnerHorizons;
    packingParam.permutations = options.usePermutations;
    packingParam.cellSize = 1;
    packingParam.rotationNum = 16; //number of rasterizations in 360°
    //packingParam.rotationTicksHalfPi = 4; //number of rasterizations in 90°

    // TODO padding computation should depend on the destination size and should be handled better, really...
    std::vector<double> containerRatios = ComputeContainerRatios(graph->textureObject, options.oneContainer);

    int sizeUnit;
    if (options.lowResPacking) {
        sizeUnit = 256;
        packingParam.pad = 2;
    } else {
        sizeUnit = 2048;
        packingParam.pad = 4;
    }

    std::vector<Point2i> containerVec;
    for (auto ratio : containerRatios) {
        Point2i gridSize(sizeUnit, sizeUnit);
        if (ratio > 1)
            gridSize.X() *= ratio;
        else if (ratio < 1)
            gridSize.Y() /= ratio;
        containerVec.push_back(gridSize);
    }

    std::vector<Similarity2f> packingTransforms;
    std::vector<int> containerIndices;

    RasterizationBasedPacker::PackingStats stats;
    RasterizationBasedPacker::Pack(texOutlines, containerVec, packingTransforms, containerIndices, packingParam, stats);

    LOG_INFO << "Packing took " << timer.TimeSinceLastCheck() << " seconds";

    for (auto entry : outlineMap) {
        for (auto fptr : graph->charts[entry.first]->fpVec) {
            int ic = containerIndices[entry.second];
            Point2i gridSize = containerVec[ic];
            for (int j = 0; j < fptr->VN(); ++j) {
                Point2d uv = fptr->WT(j).P();
                Point2f p = packingTransforms[entry.second] * (Point2f(uv[0], uv[1]));
                p.X() /= (double) gridSize.X();
                p.Y() /= (double) gridSize.Y();
                fptr->WT(j).P() = Point2d(p.X(), p.Y());
                fptr->WT(j).N() = ic;
            }
        }
    }

    if (options.lowResPacking) {
        // transform the outlines to match the computed atlas
        for (std::size_t i = 0; i < texOutlines.size(); ++i) {
            int ic = containerIndices[i];
            Point2i gridSize = containerVec[ic];
            for (std::size_t j = 0; j < texOutlines[i].size(); ++j) {
                Point2f p = packingTransforms[i] * texOutlines[i][j];
                //texOutlines[i][j] = Point2f(p.X() / gridSize.X(), p.Y() / gridSize.Y());
                int minDim = std::min(gridSize.X(), gridSize.Y());
                texOutlines[i][j] = Point2f(p.X() / minDim, p.Y() / minDim);
            }
        }

        std::vector<Similarity2f> shiftingTransforms(texOutlines.size());
        std::vector<Point2i> scaleFactors(texOutlines.size());

        // for each container, pick all the outlines mapped to that container and
        // invoke the shifting procedure
        for (unsigned c = 0; c < containerVec.size(); ++c) {
            std::vector<int> outlineIndex;
            std::vector<Outline2f> containerOutlines;
            for (unsigned i = 0; i < texOutlines.size(); ++i) {
                if (containerIndices[i] == (int)c) {
                    outlineIndex.push_back(i);
                    containerOutlines.push_back(texOutlines[i]);
                }
            }
            std::vector<Similarity2f> containerTransforms;
            //PixelShiftingOptimizer<float, QtOutline2Rasterizer>::Parameters par = {4096, 4096, 16};

            double ratio = containerRatios[c];
            int cw = 4096;
            int ch = 4096;
            if (ratio > 1)
                cw = cw * ratio;
            else if (ratio < 1)
                ch = ch / ratio;

            PixelShiftingOptimizer<float, QtOutline2Rasterizer>::Parameters par = {cw, ch, 16};


            Point2i scale = PixelShiftingOptimizer<float, QtOutline2Rasterizer>::Apply(containerOutlines, containerTransforms, par);

            for (unsigned j = 0; j < containerTransforms.size(); ++j) {
                shiftingTransforms[outlineIndex[j]] = containerTransforms[j];
                scaleFactors[outlineIndex[j]] = scale;
            }
        }

        // apply the transforms
        for (auto entry : outlineMap) {
            for (auto fptr : graph->charts[entry.first]->fpVec) {
                for (int j = 0; j < fptr->VN(); ++j) {
                    Point2d uv = fptr->WT(j).P();

                    int ic = containerIndices[entry.second];
                    Point2i gridSize = containerVec[ic];
                    int minDim = std::min(gridSize.X(), gridSize.Y());

                    uv.X() *= (gridSize.X() / (double) minDim);
                    uv.Y() *= (gridSize.Y() / (double) minDim);

                    Point2f p(uv[0], uv[1]);

                    p = shiftingTransforms[entry.second] * p;
                    p.X() /= (double) scaleFactors[entry.second].X();
                    p.Y() /= (double) scaleFactors[entry.second].Y();

                    fptr->WT(j).P() = Point2d(p.X(), p.Y());
                }
            }
        }
    }

    for (auto& entry : graph->charts)
        entry.second->ParameterizationChanged(); // packing changes texture coords even for charts that are not reparameterized

    stats.efficiency = graph->AreaUV();
    return stats;
}

