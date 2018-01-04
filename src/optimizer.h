#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/poisson_solver.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <vcg/space/rasterized_outline2_packer.h>
#include <vcg/space/outline2_packer.h>
#include <wrap/qt/outline2_rasterizer.h>

#include <utility>
#include <vector>

#include "mesh_graph.h"
#include "timer.h"
#include "dcpsolver.h"
#include "fixed_border_bijective.h"
#include <ext/texcoord_optimization.h>
#include "iterative.h"

#include "uv.h"

#include <vcg/space/index/grid_util2d.h>
#include <vcg/space/segment2.h>
#include <vcg/space/intersection2.h>


struct Point2iHasher {
    std::size_t operator()(const Point2i& p) const noexcept
    {
        std::size_t seed = 0;
        seed ^= std::hash<int>()(p[0]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<int>()(p[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

template <typename ScalarType = float>
static void MeshOutlinesUV(Mesh& m, std::vector<std::vector<Point2<ScalarType>>> &outline2Vec)
{
    tri::UpdateFlags<Mesh>::FaceClearV(m);
    tri::UpdateFlags<Mesh>::VertexClearV(m);
    tri::UpdateTopology<Mesh>::FaceFace(m);

    outline2Vec.clear();
    std::vector<Point2<ScalarType>> outline;

    for (auto& f : m.face) {
        auto fptr = &f;
        for (int i = 0; i < 3; ++i) {
            if (!fptr->IsV() && face::IsBorder(*fptr, i)) {
                face::Pos<Mesh::FaceType> p(fptr, i);
                face::Pos<Mesh::FaceType> startPos = p;
                assert(p.IsBorder());
                do {
                    assert(p.IsManifold());
                    p.F()->SetV();
                    //outline.push_back(Point2<ScalarType>(p.V()->P()));
                    Point2d uv = p.F()->WT(p.VInd()).P();
                    outline.push_back(Point2<ScalarType>{ScalarType(uv[0]), ScalarType(uv[1])});
                    p.NextB();
                }
                while (p != startPos);
                outline2Vec.push_back(outline);
                outline.clear();
            }
        }
    }
}

// Computes the texture outlines of a given chart
/// FIXME if the chart is an aggregate that has not been reparameterized this breaks...
template <typename ScalarType = float>
static void ChartOutlinesUV(Mesh& m, GraphManager::ChartHandle chart, std::vector<std::vector<Point2<ScalarType>>> &outline2Vec)
{
    struct FaceAdjacency {
        bool changed[3] = {false, false, false};
        Mesh::FacePointer FFp[3];
        char FFi[3];
    };

    auto CCIDh = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(m, CCIDh));

    std::unordered_map<Mesh::FacePointer,FaceAdjacency> savedAdj;
    for (auto fptr : chart->fpVec) {
        fptr->ClearV();
        for (int i = 0; i < 3; ++i) {
            if (CCIDh[fptr] != CCIDh[fptr->FFp(i)]) {
                savedAdj[fptr].changed[i] = true;
                savedAdj[fptr].FFp[i] = fptr->FFp(i);
                savedAdj[fptr].FFi[i] = fptr->FFi(i);
                fptr->FFp(i) = fptr;
                fptr->FFi(i) = i;
            }
        }
    }

    outline2Vec.clear();
    std::vector<Point2<ScalarType>> outline;

    for (auto fptr : chart->fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (!fptr->IsV() && face::IsBorder(*fptr, i)) {
                face::Pos<Mesh::FaceType> p(fptr, i);
                face::Pos<Mesh::FaceType> startPos = p;
                assert(p.IsBorder());
                do {
                    assert(p.IsManifold());
                    p.F()->SetV();
                    //outline.push_back(Point2<ScalarType>(p.V()->P()));
                    Point2d uv = p.F()->WT(p.VInd()).P();
                    outline.push_back(Point2<ScalarType>{ScalarType(uv[0]), ScalarType(uv[1])});
                    p.NextB();
                }
                while (p != startPos);
                outline2Vec.push_back(outline);
                outline.clear();
            }
        }
    }

    for (auto& entry : savedAdj) {
        for (int i = 0; i < 3; ++i) {
            if (entry.second.changed[i]) {
                entry.first->FFp(i) = entry.second.FFp[i];
                entry.first->FFi(i) = entry.second.FFi[i];
            }
        }
    }
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

static bool ChartParameterizationHasOverlaps(Mesh& m, GraphManager::ChartHandle chart)
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

/// TODO REFACTOR
enum DirectParameterizer {
    DCP, FixedBorderBijective
};

enum TexCoordOptimizer {
    AreaPreserving, SymmetricDirichletOpt, MIPS
};

struct ParameterizationStrategy {
    DirectParameterizer directParameterizer = DCP;
    TexCoordOptimizer optimizer = AreaPreserving;
    ParameterizationGeometry geometry = Model;
    int optimizerIterations = 0;
};

static bool ParameterizeMesh(Mesh& m, ParameterizationStrategy strategy)
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

static bool ParameterizeChartFromInitialTexCoord(Mesh &m, GraphManager::ChartHandle ch, ParameterizationStrategy strategy)
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

/*
bool ParameterizeChartFromInitialTexCoord(Mesh &m, GraphManager::ChartHandle ch)
{
    auto CCIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(m, CCIDh));

    auto ICCh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "InitialConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(m, ICCh));

    PMesh pm;
    // Count the required vertices, split those at a seam in the initial parameterization
    std::unordered_map<Mesh::VertexPointer,std::unordered_map<RegionID,PMesh::VertexPointer>> mv_to_pmv{ch->FN() * 3};

    std::size_t vn = 0;
    for (auto fptr : ch->fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (mv_to_pmv.count(fptr->V(i)) == 0 || mv_to_pmv[fptr->V(i)].count(ICCh[fptr]) == 0) {
                vn++;
                mv_to_pmv[fptr->V(i)].insert(std::make_pair(ICCh[fptr], nullptr));
            }
        }
    }

    auto pmvi = tri::Allocator<PMesh>::AddVertices(pm, vn);
    auto pmfi = tri::Allocator<PMesh>::AddFaces(pm, ch->FN());

    for (auto fptr : ch->fpVec) {
        PMesh::FacePointer pmf = &*pmfi++;
        for (int i = 0; i < 3; ++i) {
            Mesh::VertexPointer mv = fptr->V(i);
            PMesh::VertexPointer& pmv = mv_to_pmv[mv][ICCh[fptr]];
            if (pmv == nullptr) pmv = &*pmvi++;
            pmf->V(i) = pmv;

            // TODO maybe add a parameter to choose the parameterization source (mesh geometry or old parameterization)
            pmv->P() = PMesh::VertexType::CoordType{fptr->cWT(i).U(), fptr->cWT(i).V(), 0};
            //pmv->P() = fptr->V(i)->P();
        }
    }

    for (auto fptr : ch->fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (face::IsBorder(*fptr, i) || CCIDh[fptr] != CCIDh[fptr->FFp(i)]) { // Border in texture space
                for (auto& idVertPair : mv_to_pmv[fptr->V0(i)]) idVertPair.second->SetB();
                for (auto& idVertPair : mv_to_pmv[fptr->V1(i)]) idVertPair.second->SetB();
            }
        }
    }

    DCPSolver<PMesh> solver(pm);
    for (auto& entry : mv_to_pmv) {
        if (entry.second.size() > 1) {
            std::vector<PMesh::VertexPointer> vertGroup;
            for (auto& idVertPair : entry.second) {
                vertGroup.push_back(idVertPair.second);
            }
            solver.DeclareUnique(vertGroup);
        }
    }

    // TEST CODE
    bool solved = solver.Solve(DCPSolver<PMesh>::__trust_vertex_border_flags);
    //bool solved = solver.Solve();

    if (solved) { // Copy texture coords back
        for (auto& fptr : ch->fpVec) {
            for (int i = 0; i < 3; ++i) {
                auto vi = fptr->V(i);
                assert(mv_to_pmv[vi].size() > 0);
                auto pmv = mv_to_pmv[vi].begin()->second;
                fptr->WT(i).P() = pmv->T().P();
            }

        }
    }

#ifndef NDEBUG

    if (!solved) std::cout << "Not solved" << std::endl;

    int n1 = 0;
    int n2 = 0;
    int n3 = 0;

    for (auto &f : pm.face) {
        Triangle3<float> tri {
            Point3f { f.WT(0).P().X(), f.WT(0).P().Y(), 0.0f },
            Point3f { f.WT(1).P().X(), f.WT(1).P().Y(), 0.0f },
            Point3f { f.WT(2).P().X(), f.WT(2).P().Y(), 0.0f }
        };
        float area = vcg::DoubleArea(tri);
        if (area < 0.0f) n1++;
        else if (area > 0.0f) n2++;
        else n3++;
    }

    std::cout << "Negative area " << n1 << std::endl;
    std::cout << "Positive area " << n2 << std::endl;
    std::cout << "Zero area     " << n3 << std::endl;

    //tri::io::ExporterOBJ<PMesh>::Save(pm, "test_mesh.obj", tri::io::Mask::IOM_WEDGTEXCOORD);

#endif

    return solved;
}
*/


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
static int ParameterizeGraph(GraphManager& gm, ParameterizationStrategy strategy = ParameterizationStrategy{},
                             bool failsafe = true)
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






#if 0
// returns the number of charts that could not be parameterized
/// TODO update distortion info if needed (this should also be done through the graph manager)
static int ParameterizeGraph(std::shared_ptr<MeshGraph> graph, ParameterizationStrategy strategy = ParameterizationStrategy{})
{
    Timer timer;
    Mesh& m = graph->mesh;

    std::cout << "Parameterizing " << graph->charts.size() << " regions..." << std::endl;

    std::vector<std::vector<Point2f>> texOutlines;
    std::unordered_map<RegionID,std::size_t> outlineMap; // map each region to the index of its outline in texOutlines
    int iter = 0;
    int numFailed = 0;
    for (auto entry : graph->charts) {
        std::shared_ptr<FaceGroup> chart = entry.second;

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

            if (ChartParameterizationHasOverlaps(m, chart)) {
                std::cout << "WARNING: CHART " << chart->id << " HAS OVERLAPS IN THE PARAMETERIZATION" << std::endl;
                // store this chart in a list and split it later to avoid invalidating the iterator
            }

            if (parameterized) {
                // Normalize area: the region gets scaled by sqrt(oldUVArea)/sqrt(newUVArea) to keep the original proportions
                // between the regions of the atlas in an attempt to not lose too much texture data when the new texture is rendered
                float newUvArea = chart->AreaUV();
                float scale = std::sqrt(oldUvArea / newUvArea);
                assert(scale > 0);
                vcg::Box2d uvBox = chart->UVBox();
                for (auto fptr : chart->fpVec) {
                    for (int i = 0; i < 3; ++i) {
                        fptr->WT(i).P() = (fptr->WT(i).P() - uvBox.min) * scale;
                    }
                }

                chart->numMerges = 0;
            }
            else {
                numFailed++;
                std::cout << "Parameterization failed" << std::endl;
            }
        }
        // Save the outline of the parameterization for this portion of the mesh

        std::vector<std::vector<Point2f>> uvOutlines;
        ChartOutlinesUV(m, chart, uvOutlines);
        int i = tri::OutlineUtil<float>::LargestOutline2(uvOutlines);
        if (tri::OutlineUtil<float>::Outline2Area(uvOutlines[i]) < 0)
            tri::OutlineUtil<float>::ReverseOutline2(uvOutlines[i]);

        outlineMap[chart->id] = texOutlines.size();
        texOutlines.push_back(uvOutlines[i]);

        std::cout << "Iteration " << iter++ << " took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
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
#endif

/// TODO refactoring
static void ReduceTextureFragmentation_NoPacking(GraphManager &gm, std::size_t minRegionSize);

static void ReduceTextureFragmentation(Mesh &m, std::shared_ptr<MeshGraph> graph, std::size_t minRegionSize)
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

static void ReduceTextureFragmentation_NoPacking(GraphManager &gm, std::size_t minRegionSize)
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

#if 0
static void ReduceTextureFragmentation2(Mesh &m, std::shared_ptr<MeshGraph> graph, std::size_t minRegionSize)
{
    if (minRegionSize == 0) return;

    assert(minRegionSize < (std::size_t) m.FN());

    tri::UpdateTopology<Mesh>::FaceFace(m);

    GraphManager gm{graph};

    std::cout << "Initialized graph manager" << std::endl;

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

    timer.Reset();

    std::cout << "Parameterizing " << graph->charts.size() << " regions..." << std::endl;
    int iter = 0;

    std::vector<std::vector<Point2f>> texOutlines;
    std::unordered_map<std::size_t, std::size_t> outlineMap;
    // reparameterize each region with the poisson solver
    for (auto entry : graph->charts) {
        std::shared_ptr<FaceGroup> chart = entry.second;

        std::cout << "Chart " << chart->id << " - FN=" << chart->FN() << ", FI=" << tri::Index(m, chart->Fp()) << std::endl;

        // If the old chart has no valid parameterization, skip it
        float oldUvArea = chart->AreaUV();
        if (oldUvArea == 0) continue;

        // Build the mesh to parameterize
        Timer t_;
        PMesh pm;
        for (auto fptr : chart->fpVec) {
            tri::Allocator<PMesh>::AddFace(pm, fptr->P(0), fptr->P(1), fptr->P(2));
        }
        tri::Clean<PMesh>::RemoveDuplicateVertex(pm);
        tri::Allocator<PMesh>::CompactEveryVector(pm);
        std::cout << "Per face mesh alloc took " << t_.TimeElapsed() << std::endl;
        //tri::UpdateTopology<PMesh>::FaceFace(pm);
        tri::UpdateSelection<PMesh>::Clear(pm);

        //tri::io::ExporterOBJ<PMesh>::Save(pm, "submesh_iteration.obj", tri::io::Mask::IOM_WEDGTEXCOORD);

        tri::PoissonSolver<PMesh> solver(pm);
        //DCPSolver<PMesh> solver(pm);
        tri::UpdateBounding<PMesh>::Box(pm);
        if (!solver.IsFeasible()) {

            // Note: since we merge two regions only if the result is parameterizable
            // with the PoissonSolver, unfeasibility can only occur for regions that have NOT been merged with anything
            // else, thus having the original parameterization untouched. Just assert the condition and copy the original
            // parameterization over to the PMesh

            assert(chart->numMerges == 0);

            for (std::size_t i = 0; i < pm.face.size(); ++i) {
                auto &pf = pm.face[i];
                auto &f = *(chart->fpVec[i]);
                for (int k = 0; k < f.VN(); ++k) {
                    pf.WT(k) = f.WT(k);
                }
            }

            // !!! this does not happen to work as there are issues with selecting the outline2 from such regions
            // Example: closed surfaces parameterized without cuts (like two sides projected to a plane)

        } else {

            // In this case generate a new parameterization

            // Solve the system and update texture coords
            solver.Init();
            solver.FixDefaultVertices();
            solver.SolvePoisson();
            //tri::UpdateTexture<PMesh>::WedgeTexFromVertexTex(pm);

            // Normalize area: the region gets scaled by sqrt(oldUVArea)/sqrt(newUVArea) to keep the original proportions
            // between the regions of the atlas in an attempt to not lose too much texture data when the new texture is rendered
            float newUvArea = 0.0f;
            for (auto &pf : pm.face) {
                newUvArea += std::abs(tri::Distortion<PMesh,true>::AreaUV(&pf));
            }

            float scale = std::sqrt(oldUvArea / newUvArea);

            /*if (scale <= 0) {
                char name[50] = {0};
                std::sprintf(name, "chart_%d.obj", iter);
                tri::io::ExporterOBJ<PMesh>::Save(pm, name, tri::io::Mask::IOM_WEDGTEXCOORD);
            }*/

            assert(scale > 0);

            vcg::Box2f uvBox = tri::UV_Utils<PMesh>::PerWedgeUVBox(pm);
            for (std::size_t i = 0; i < pm.face.size(); ++i) {
                auto &pf = pm.face[i];
                auto &f = *(chart->fpVec[i]);
                for (int k = 0; k < f.VN(); ++k) {
                    pf.WT(k).P() = (pf.WT(k).P() - uvBox.min) * scale;
                    f.WT(k) = pf.WT(k);
                }
            }
        }

        // Save the outline of the parameterization for this portion of the mesh

        std::vector<std::vector<Point2f>> regionOutlines;
        ConvertTextureBoundaryToOutline2Vec<float, PMesh>(pm, regionOutlines);
        int i = tri::OutlineUtil<float>::LargestOutline2(regionOutlines);
        if (tri::OutlineUtil<float>::Outline2Area(regionOutlines[i]) < 0)
            tri::OutlineUtil<float>::ReverseOutline2(regionOutlines[i]);

        outlineMap[entry.first] = texOutlines.size();
        texOutlines.push_back(regionOutlines[i]);

        std::cout << "Iteration " << iter++ << " took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    }

    std::cout << "Parameterization took " << timer.TimeElapsed() << " seconds" << std::endl;

    std::cout << "Packing the atlas..." << std::endl;

    // pack the atlas
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
                Point2f transformedTexCoordPos = transforms[p.second] * (fptr->WT(j).P());
                fptr->WT(j).P() = transformedTexCoordPos / 1024.0f;
            }
        }
    }

}
#endif

#endif // OPTIMIZER_H

