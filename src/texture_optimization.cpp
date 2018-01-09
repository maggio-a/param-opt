#include "texture_optimization.h"

#include "uv.h"
#include "mesh_graph.h"
#include "timer.h"
#include "dcpsolver.h"
#include "fixed_border_bijective.h"
#include <ext/texcoord_optimization.h>
#include "iterative.h"

#include <vcg/complex/complex.h>
#include <vcg/space/rasterized_outline2_packer.h>
#include <vcg/space/outline2_packer.h>
#include <wrap/qt/outline2_rasterizer.h>
#include <vcg/space/index/grid_util2d.h>
#include <vcg/space/segment2.h>
#include <vcg/space/intersection2.h>

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

bool ParameterizeMesh(Mesh& m, ParameterizationStrategy strategy)
{
    WedgeTexCoordAttributePosition<Mesh> texcoordPosition{m, "WedgeTexCoordStorage"};
    bool solved = false;
    if (strategy.directParameterizer == DirectParameterizer::DCP) {
        DCPSolver<Mesh> solver(m);
        if (strategy.geometry == ParameterizationGeometry::Model)
            solved = solver.Solve(DefaultVertexPosition<Mesh>{});
        else
            solved = solver.Solve(texcoordPosition);
    } else if (strategy.directParameterizer == DirectParameterizer::FixedBorderBijective) {
        UniformSolver<Mesh> solver(m);
        solved = solver.Solve();
    } else {
        assert(0 && "Solver not supported");
    }

    if (solved) { // optimize and Copy texture coords back

        if (strategy.optimizerIterations > 0) {

            //GradientDescent opt(std::make_shared<SymmetricDirichlet>(pm));
            Energy::Geometry mode = Energy::Geometry::Model;
            if (strategy.geometry == ParameterizationGeometry::Texture) mode = Energy::Geometry::Texture;

            LBFGS opt(std::make_shared<SymmetricDirichlet>(m, mode), 10);
            //GradientDescent opt(std::make_shared<SymmetricDirichlet>(m, mode));

            Timer t;
            if (strategy.optimizerIterations < 0) {
                std::abort();
            } else {
                int i;
                double energyVal, gradientNorm, energyDiff;
                for (i = 0; i < strategy.optimizerIterations; ++i) {
                    energyVal = opt.Iterate(gradientNorm, energyDiff);
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
            }
            std::cout << "Optimization took " << t.TimeSinceLastCheck() << " seconds" << std::endl;
        } else {
          /*  strategy.optimizerIterations = -strategy.optimizerIterations;

            tri::TexCoordOptimization<Mesh> *opt;
            opt = new tri::MIPSTexCoordOptimization<Mesh, WedgeTexCoordAttributePosition<Mesh>>(m, texcoordPosition);

            opt->SetNothingAsFixed();
            opt->TargetCurrentGeometry();

            Timer t;
            if (strategy.optimizerIterations < 0) {
            } else {
                float movement;
                for (int i = 0; i < strategy.optimizerIterations; ++i) {
                    movement = opt->Iterate();
                }
                std::cout << "Max movement at last step was " << movement << std::endl;
            }
            std::cout << "Optimization took " << t.TimeSinceLastCheck() << " seconds" << std::endl;
            delete opt;*/
        }
    }
    return solved;
}

bool ParameterizeChartFromInitialTexCoord(Mesh &m, GraphManager::ChartHandle ch, ParameterizationStrategy strategy)
{
    Mesh pm;
    std::unordered_map<Mesh::VertexPointer,Mesh::VertexPointer> mv_to_pmv;

    CopyFaceGroupIntoMesh(pm, *ch, mv_to_pmv);

    bool solved = ParameterizeMesh(pm, strategy);

    if (solved) {
        for (auto& fptr : ch->fpVec) {
            for (int i = 0; i < 3; ++i) {
                auto vi = fptr->V(i);
                auto pmv = mv_to_pmv[vi];
                fptr->WT(i).P() = pmv->T().P();
            }
        }
    }

    return solved;
}


// returns the number of charts that could not be parameterized
/// TODO update distortion info if needed (this should also be done through the graph manager)
/*
 * Parameterize the mesh graph. Each region is a connected set of mesh faces, and it is assumed to be homeomorphic to a disk.
 * The parameterization of each region is performed according to the parameterization strategy passed as parameter.
 * If failsafe is true, then each region is tested for overlaps in the parameterization: in case the parameterization contains
 * overlaps, the region is split in its original components and each face is assigned its original texture coordinates. This
 * ensures that no overlaps are introduced by this procedure, potentially reducing the whole procedure to a no-op if necessary
 * (that is, if every region parameterization contains overlaps).
 *
 * After each region is parameterized the procedure packs the texture atlas.
 *
 */
int ParameterizeGraph(GraphManager& gm,
                      ParameterizationStrategy strategy,
                      bool failsafe)
{
    Timer timer;
    auto graph = gm.Graph();
    Mesh& m = graph->mesh;

    std::cout << "Parameterizing " << graph->charts.size() << " regions..." << std::endl;

    std::set<GraphManager::ChartHandle> toSplit;
    int iter = 0;
    int numFailed = 0;
    for (auto entry : graph->charts) {
        GraphManager::ChartHandle chart = entry.second;

        chart->ParameterizationChanged(); // always changes due to repacking, even if the original coordinates were preserved

        std::cout << "Chart " << chart->id << " - FN=" << chart->FN() << ", FI=" << tri::Index(m, chart->Fp()) << std::endl;

        // If the old chart has no valid parameterization, skip it
        float oldUvArea = chart->AreaUV();
        if (oldUvArea == 0) continue;

        if (chart->numMerges == 0) {
            // The chart was not merged, therefore there is nothing to do -- the original parameterization is preserved

            // !!! this does not happen to work as there are issues with selecting the outline2 from such regions
            // Example: closed surfaces parameterized without cuts (like two sides projected to a plane)

            // TODO a possible workaround would be to use the bounding box as the outline
        }
        else {
            bool parameterized = ParameterizeChartFromInitialTexCoord(m, chart, strategy);

            if (failsafe) {
                if (ChartParameterizationHasOverlaps(m, chart)) {
                    // store this chart in a list and split it later to avoid invalidating the iterator
                    std::cout << "WARNING: REGION " << chart->id << " HAS OVERLAPS IN THE PARAMETERIZATION" << std::endl;
                    toSplit.insert(chart);
                }
            }

            if (parameterized) {
                // if the parameterization is valid, scale chart relative to the original parameterization area value
                if (toSplit.find(chart) == toSplit.end()) {
                    // Normalize area: the region gets scaled by sqrt(oldUVArea)/sqrt(newUVArea) to keep the original proportions
                    // between the regions of the atlas in an attempt to not lose too much texture data when the new texture is rendered

                    float newUvArea = chart->AreaUV();
                    float scale = std::sqrt(oldUvArea / newUvArea);
                    assert(scale > 0);
                    std::cout << "Scale value = " << scale << std::endl; // shoud be very close to 1 if we optimize for area distortion
                    vcg::Box2d uvBox = chart->UVBox();
                    for (auto fptr : chart->fpVec) {
                        for (int i = 0; i < 3; ++i) {
                            fptr->WT(i).P() = (fptr->WT(i).P() - uvBox.min) * scale;
                        }
                    }

                    chart->numMerges = 0;
                }
            }
            else {
                numFailed++;
                std::cout << "Parameterization failed" << std::endl;
            }
        }

        std::cout << "Iteration " << iter++ << " took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    }

    // Handle chart that need to be split, if any
    auto WTCSh = tri::Allocator<Mesh>::FindPerFaceAttribute<TexCoordStorage>(m, "WedgeTexCoordStorage");
    assert(tri::Allocator<Mesh>::IsValidHandle<TexCoordStorage>(m, WTCSh));

    for (auto chart : toSplit) {
        // restore original texture coordinates for each face of the chart
        std::cout << "Splitting region " << chart->id << std::endl;
        for (auto fptr : chart->fpVec) {
            TexCoordStorage tcs = WTCSh[fptr];
            for (int i = 0; i < 3; ++i) {
                fptr->WT(i) = tcs.tc[i];
                fptr->V(i)->T() = tcs.tc[i];
            }
        }

        // split the chart in the mesh graph using the graph manager
        gm.Split(chart->id);
    }

    // Pack the atlas
    std::vector<std::vector<Point2f>> texOutlines;
    std::unordered_map<RegionID,std::size_t> outlineMap; // map each region to the index of its outline in texOutlines

    for (auto entry : graph->charts) {
        GraphManager::ChartHandle chart = entry.second;

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

    //RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Pack(texOutlines, gridSize, transforms, packingParam);

    Point2f cover;
    //PolyPacker<float>::PackAsObjectOrientedRect(texOutlines, gridSize, transforms, cover);
    PolyPacker<float>::PackAsAxisAlignedRect(texOutlines, gridSize, transforms, cover);

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

    int c = ParameterizeGraph(gm);
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