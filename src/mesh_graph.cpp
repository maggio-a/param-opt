#include "mesh_graph.h"

#include "mesh.h"
#include "gl_utils.h"
#include "math_utils.h"
#include "metric.h"
#include "uv.h"
#include "mesh_attribute.h"
#include "uniform_solver.h"
#include "mean_value_param.h"
#include "timer.h"

#include "polygon2_triangulator.h"

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>
#include <vcg/math/histogram.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/attribute_seam.h>

#include <wrap/io_trimesh/export.h>

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>

#include <QImage>

static constexpr int MAX_REMESH_ITER = 6;

inline Polyline2 BuildPolyline2(const std::vector<std::size_t> &vfi, const std::vector<int> &vvi, const Mesh& shell2D)
{
    Polyline2 polyline;
    polyline.reserve(vfi.size());
    for (std::size_t i = 0; i < vfi.size(); ++i) {
        const vcg::Point3d& p = shell2D.face[vfi[i]].cP(vvi[i]);
        assert(p.Z() == 0.0);
        polyline.push_back({p.X(), p.Y()});
    }
    return polyline;
}

// assumes topology is correct
void ComputeBoundaryInfo(Mesh& m)
{
    BoundaryInfo& info = GetBoundaryInfoAttribute(m)();
    info.Clear();
    tri::UpdateFlags<Mesh>::FaceClearV(m);
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            if (!f.IsV() && face::IsBorder(f, i)) {
                double totalBorderLength = 0;
                std::vector<std::size_t> borderFaces;
                std::vector<int> vi;

                face::Pos<Mesh::FaceType> p(&f, i);
                face::Pos<Mesh::FaceType> startPos = p;
                assert(p.IsBorder());
                do {
                    assert(p.IsManifold());
                    p.F()->SetV();
                    borderFaces.push_back(tri::Index(m, p.F()));
                    vi.push_back(p.VInd());
                    totalBorderLength += EdgeLength(*p.F(), p.VInd());
                    p.NextB();
                } while (p != startPos);
                info.vBoundaryLength.push_back(totalBorderLength);
                info.vBoundarySize.push_back(borderFaces.size());
                info.vBoundaryFaces.push_back(borderFaces);
                info.vVi.push_back(vi);
            }
        }
    }

    std::cout << "Mesh has " << info.N() << " boundaries" << std::endl;
}

void CloseMeshHoles(Mesh& shell)
{
    int startFN = shell.FN();

    // Get border info
    assert(HasBoundaryInfoAttribute(shell));
    BoundaryInfo& info = GetBoundaryInfoAttribute(shell)();

    // Leave only the longest boundary
    tri::UpdateFlags<Mesh>::FaceClearS(shell);
    assert(info.vBoundaryFaces.size() > 0 && "Mesh has no boundaries");
    if (info.vBoundaryFaces.size() > 1) {
        std::size_t k = info.LongestBoundary();
        // select all the boundary faces
        for (std::size_t i = 0; i < info.vBoundaryFaces.size(); ++i) {
            if (i == k) continue;
            for (auto j : info.vBoundaryFaces[i]) {
                assert(face::IsBorder(shell.face[j], 0) || face::IsBorder(shell.face[j], 1) || face::IsBorder(shell.face[j], 2));
                shell.face[j].SetS();
            }
        }
        tri::Hole<Mesh>::EarCuttingFill<tri::MinimumWeightEar<Mesh>>(shell, shell.FN(), true);
    }
    tri::UpdateTopology<Mesh>::FaceFace(shell);

    auto ia = GetFaceIndexAttribute(shell);
    tri::UpdateFlags<Mesh>::FaceClearS(shell);
    for (auto& f: shell.face) {
        if (int(tri::Index(shell, f)) >= startFN) {
            f.SetHoleFilling();
            ia[f] = -1;
        }
    }
}

void ColorFace(Mesh& shell)
{
    for (auto& sf : shell.face) {
        if (!sf.IsMesh()) {
            sf.C() = vcg::Color4b::Cyan;
        } else {
            sf.C() = vcg::Color4b::LightGreen;
        }
    }
}

void PrintParameterizationInfo(std::shared_ptr<MeshGraph> pdata)
{
    std::cout << "Parameterization Info (Charts, Area 3D, Area UV, Border UV) " << std::endl;
    std::cout << pdata->charts.size() << " " << pdata->Area3D() << " "
              << pdata->AreaUV() << " " << pdata->BorderUV() << std::endl;
}

std::shared_ptr<MeshGraph> ComputeParameterizationGraph(Mesh &m, TextureObjectHandle textureObject)
{
    std::size_t numRegions = ComputePerFaceConnectedComponentIdAttribute(m);

    std::shared_ptr<MeshGraph> graph = std::make_shared<MeshGraph>(m);
    graph->textureObject = textureObject;
    graph->charts.reserve(numRegions);

    auto CCIDh = GetConnectedComponentIDAttribute(m);

    auto ICCIDh  = GetInitialConnectedComponentIDAttribute(m);
    ICCIDh._handle->data.assign(CCIDh._handle->data.begin(), CCIDh._handle->data.end());

    // build parameterization graph
    tri::UpdateTopology<Mesh>::FaceFace(m);
    for (auto &f : m.face) {
        RegionID regionId = CCIDh[&f];
        graph->GetChart_Insert(regionId)->AddFace(&f);
        // TODO this may be refactored into AddFace
        for (int i = 0; i < f.VN(); ++i) {
            RegionID adjId = CCIDh[f.FFp(i)];
            if (regionId != adjId) {
                (graph->GetChart_Insert(regionId)->adj).insert(graph->GetChart_Insert(adjId));
            }
        }
    }

    return graph;
}

void CopyToMesh(FaceGroup& fg, Mesh& m)
{
    m.Clear();
    auto ia = GetFaceIndexAttribute(m);
    std::unordered_map<Mesh::VertexPointer, Mesh::VertexPointer> vpmap;
    vpmap.reserve(fg.FN() * 3);
    std::size_t vn = 0;
    for (auto fptr : fg.fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (vpmap.count(fptr->V(i)) == 0) {
                vn++;
                vpmap[fptr->V(i)] = nullptr;
            }
        }
    }
    auto mvi = tri::Allocator<Mesh>::AddVertices(m, vn);
    auto mfi = tri::Allocator<Mesh>::AddFaces(m, fg.FN());
    for (auto fptr : fg.fpVec) {
        Mesh::FacePointer mfp = &*mfi++;
        ia[mfp] = tri::Index(fg.mesh, fptr);
        for (int i = 0; i < 3; ++i) {
            Mesh::VertexPointer vp = fptr->V(i);
            typename Mesh::VertexPointer& mvp = vpmap[vp];
            if (mvp == nullptr) {
                mvp = &*mvi++;
                mvp->P() = vp->P();
            }
            mfp->V(i) = mvp;
            mfp->WT(i) = fptr->WT(i); // TODO should this be the stored initial wt?
        }
        mfp->SetMesh();
    }

    std::cout << "Built mesh has " << m.FN() << " faces" << std::endl;
}

void SyncShellWithUV(Mesh& shell)
{
    for (auto& v : shell.vert) {
        v.P().X() = v.T().U();
        v.P().Y() = v.T().V();
        v.P().Z() = 0.0;
    }
}

void SyncShellWithModel(Mesh& shell, Mesh& baseMesh)
{
    auto ia = GetFaceIndexAttribute(shell);
    for (auto& sf : shell.face) {
        assert(sf.IsMesh());
        auto& f = baseMesh.face[ia[sf]];
        for (int i = 0; i < 3; ++i)
            sf.P(i) = f.P(i);
    }
}

void ClearHoleFillingFaces(Mesh& shell, bool holefill, bool scaffold)
{
    for (auto& f : shell.face)
        if ((holefill && f.IsHoleFilling()) || (scaffold && f.IsScaffold()))
            tri::Allocator<Mesh>::DeleteFace(shell, f);

    tri::Clean<Mesh>::RemoveUnreferencedVertex(shell);
    tri::UpdateTopology<Mesh>::FaceFace(shell);
    tri::UpdateTopology<Mesh>::VertexFace(shell);
    tri::Allocator<Mesh>::CompactEveryVector(shell);
}

void CloseShellHoles(Mesh& shell, ParameterizationGeometry targetGeometry, Mesh& inputMesh)
{
    auto ia = GetFaceIndexAttribute(shell);

    // values needed to compute target shapes
    auto psi = GetParameterizationScaleInfoAttribute(inputMesh);
    auto tsa = GetTargetShapeAttribute(shell);
    double avg3D = psi().surfaceArea / psi().numNonZero;
    double avgUV = psi().parameterArea / psi().numNonZero;


    BoundaryInfo& info = GetBoundaryInfoAttribute(shell)();
    for (std::size_t bi = 0; bi < info.vBoundaryFaces.size(); ++bi) {
        if (bi == info.LongestBoundary())
            continue;

        Polyline2 outline = BuildPolyline2(info.vBoundaryFaces[bi], info.vVi[bi], shell);
        unsigned szoutline = outline.size();

        Poly2 poly = { outline };

        Polyline2 points;
        std::vector<unsigned> indices;

        Triangulate(poly, {}, points, indices);

        assert(points.size() >= szoutline);

        unsigned nv = points.size() - szoutline;
        unsigned nf = indices.size() / 3;

        auto vi = tri::Allocator<Mesh>::AddVertices(shell, nv);
        auto fi = tri::Allocator<Mesh>::AddFaces(shell, nf);

        std::vector<Mesh::VertexPointer> vp;
        vp.reserve(szoutline + nv);
        for (unsigned i = 0; i < szoutline; ++i)
            vp.push_back(shell.face[info.vBoundaryFaces[bi][i]].V(info.vVi[bi][i]));
        for (unsigned i = 0; i < nv; ++i) {
            (vi+i)->P()     = Point3d(points[szoutline + i][0], points[szoutline + i][1], 0);
            (vi+i)->T().P() = Point2d(points[szoutline + i][0], points[szoutline + i][1]);
            vp.push_back(&*(vi + i));
        }

        int k = 0;
        while (fi != shell.face.end()) {
            for (int i = 0; i < 3; ++i) {
                auto vi = vp[indices[k++]];
                fi->V(i) = vi;
                fi->WT(i) = vi->T();
            }
            ia[fi] = -1;
            fi->SetHoleFilling();

            CoordStorage target;
            double scale = 1.0;
            if (targetGeometry == Model)
                scale = std::sqrt(avg3D / DistortionMetric::Area3D(*fi));
            else if (targetGeometry == Texture)
                scale = std::sqrt(avgUV / DistortionMetric::Area3D(*fi));
            else
                assert(0 && "Unexpected targetGeometry parameter value");
            assert(scale > 0);
            target.P[0] = fi->P(0) * scale;
            target.P[1] = fi->P(1) * scale;
            target.P[2] = fi->P(2) * scale;
            tsa[fi] = target;

            fi++;
        }
    }

    tri::UpdateTopology<Mesh>::FaceFace(shell);
}

void RebuildScaffold(Mesh& shell, ParameterizationGeometry targetGeometry, Mesh& inputMesh)
{
    for (auto& sf : shell.face)
        if (sf.IsScaffold())
            tri::Allocator<Mesh>::DeleteFace(shell, sf);

    tri::Clean<Mesh>::RemoveUnreferencedVertex(shell);
    tri::UpdateTopology<Mesh>::FaceFace(shell);
    tri::UpdateTopology<Mesh>::VertexFace(shell);
    tri::Allocator<Mesh>::CompactEveryVector(shell);

    BuildScaffold(shell, targetGeometry, inputMesh);

    /*
    auto ia = GetFaceIndexAttribute(shell);
    auto tsa = GetTargetShapeAttribute(shell);
    int k = 0;
    for (auto& sf : shell.face) {
        if (sf.IsHoleFilling()) {
            ia[sf] = -1;
            // for scaffold faces the target shape is the initial shape
            CoordStorage target;
            target.P[0] = sf.P(0);
            target.P[1] = sf.P(1);
            target.P[2] = sf.P(2);
            tsa[sf] = target;
        }
    }
    */

}

void BuildScaffold(Mesh& shell, ParameterizationGeometry targetGeometry, Mesh& inputMesh)
{
    // Compute uv box and make it square and larger to fit the shell in
    Box2d uvBox = UVBoxVertex(shell);
    uvBox.MakeSquare();
    Point2d diag = uvBox.max - uvBox.min;
    Point2d newMax = uvBox.min + 1.15*diag;
    Point2d newMin = uvBox.max - 1.15*diag;
    uvBox.Add(newMax);
    uvBox.Add(newMin);

    BoundaryInfo& info = GetBoundaryInfoAttribute(shell)();

    Poly2Point c0 = { uvBox.min.X(), uvBox.min.Y() };
    Poly2Point c1 = { uvBox.max.X(), uvBox.min.Y() };
    Poly2Point c2 = { uvBox.max.X(), uvBox.max.Y() };
    Poly2Point c3 = { uvBox.min.X(), uvBox.max.Y() };

    constexpr int SUBDIVISION = 4;
    Polyline2 corners(4*SUBDIVISION);
    for (int i = 0; i < SUBDIVISION; ++i) {
        double t = i / (double) SUBDIVISION;
        corners[0 * SUBDIVISION + i] = Poly2PointLerp(c0, c1, t);
        corners[1 * SUBDIVISION + i] = Poly2PointLerp(c1, c2, t);
        corners[2 * SUBDIVISION + i] = Poly2PointLerp(c2, c3, t);
        corners[3 * SUBDIVISION + i] = Poly2PointLerp(c3, c0, t);

    }

    std::size_t li = info.LongestBoundary();
    Polyline2 outline = BuildPolyline2(info.vBoundaryFaces[li], info.vVi[li], shell);
    unsigned szoutline = outline.size();

    Poly2 scaffoldPoly = { outline, corners };

    Polyline2 points;
    std::vector<unsigned> indices;

    std::vector<double> holes;
    Point3d holeMarker = Barycenter(shell.face[0]);
    holes.push_back(holeMarker[0]);
    holes.push_back(holeMarker[1]);

    Triangulate(scaffoldPoly, holes, points, indices);

    // allocate new primitives
    assert(points.size() >= szoutline);
    unsigned nv = points.size() - szoutline;
    unsigned nf = indices.size() / 3;
    auto vi = tri::Allocator<Mesh>::AddVertices(shell, nv);
    auto fi = tri::Allocator<Mesh>::AddFaces(shell, nf);

    // vector of vertex pointers, indexed by 'indices'
    std::vector<Mesh::VertexPointer> vp;
    vp.reserve(szoutline + nv);
    for (unsigned i = 0; i < szoutline; ++i)
        vp.push_back(shell.face[info.vBoundaryFaces[li][i]].V(info.vVi[li][i]));
    for (unsigned i = 0; i < nv; ++i) {
        (vi+i)->P()     = Point3d(points[szoutline + i][0], points[szoutline + i][1], 0);
        (vi+i)->T().P() = Point2d(points[szoutline + i][0], points[szoutline + i][1]);
        vp.push_back(&*(vi + i));
    }

    auto ia = GetFaceIndexAttribute(shell);
    auto tsa = GetTargetShapeAttribute(shell);
    int k = 0;
    while (fi != shell.face.end()) {
        // fill face attributes
        for (int i = 0; i < 3; ++i) {
            auto vi = vp[indices[k++]];
            fi->V(i) = vi;
            fi->WT(i) = vi->T();
        }
        ia[fi] = -1;
        fi->SetScaffold();

        // for scaffold faces the target shape is the initial shape
        CoordStorage target;
        target.P[0] = fi->P(0);
        target.P[1] = fi->P(1);
        target.P[2] = fi->P(2);
        tsa[fi] = target;

        fi++;
    }

    tri::UpdateTopology<Mesh>::FaceFace(shell);
    tri::UpdateTopology<Mesh>::VertexFace(shell);
}

bool BuildShell(Mesh& shell, FaceGroup& fg, ParameterizationGeometry targetGeometry, bool useExistingUV)
{
    Mesh& m = fg.mesh;

    CopyToMesh(fg, shell);

    tri::UpdateTopology<Mesh>::FaceFace(shell);
    int splitCount = tri::Clean<Mesh>::SplitNonManifoldVertex(shell, 0.15);
    if (splitCount > 0) {
        std::cout << "Mesh was not vertex-manifold, " << splitCount << " vertices split" << std::endl;
    }
    auto ia = GetFaceIndexAttribute(shell);

    if (useExistingUV) {
        // copy wedge texture data from the input mesh
        for (auto& sf : shell.face) {
            auto& f = m.face[ia[sf]];
            for (int i = 0; i < 3; ++i) {
                sf.WT(i).P() = f.WT(i).P();
            }
        }


        // split vertices at seams
        int vn = shell.VN();
        auto vExt = [](const Mesh& msrc, const MeshFace& f, int k, const Mesh& mdst, MeshVertex& v) {
            (void) msrc;
            (void) mdst;
            v.ImportData(*(f.cV(k)));
            v.T() = f.cWT(k);
        };
        auto vCmp = [](const Mesh& mdst, const MeshVertex& v1, const MeshVertex& v2) {
            (void) mdst;
            return v1.T() == v2.T();
        };
        tri::AttributeSeam::SplitVertex(shell, vExt, vCmp);
        if (shell.VN() != (int) shell.vert.size()) {
            tri::Allocator<Mesh>::CompactEveryVector(shell);
            tri::UpdateTopology<Mesh>::FaceFace(shell);
            tri::UpdateTopology<Mesh>::VertexFace(shell);
        }

        // sync shell
        // FIXME (?) there is no guarantee that the 'outer' boundary is the same computed previously...
        for (auto& sf : shell.face) {
            for (int i = 0; i < 3; ++i) {
                sf.V(i)->T() = sf.WT(i);
            }
        }
        SyncShellWithUV(shell);
        ComputeBoundaryInfo(shell);
        CloseShellHoles(shell, targetGeometry, m);
    } else {
        // Compute the initial configuration (Tutte's parameterization)
        ComputeBoundaryInfo(shell);
        CloseMeshHoles(shell);
        UniformSolver<Mesh> solver(shell);
        bool solved = solver.Solve();
        //MeanValueSolver<Mesh> solver(shell);
        //bool solved = solver.Solve(DefaultVertexPosition<Mesh>{});
        if (!solved) {
            return false;
        } else
            SyncShellWithUV(shell);
    }

    // Compute average areas
    auto psi = GetParameterizationScaleInfoAttribute(m);
    double avg3D = psi().surfaceArea / psi().numNonZero;
    double avgUV = psi().parameterArea / psi().numNonZero;

    // Compute the target shapes
    auto tsa = GetTargetShapeAttribute(shell);
    auto wtcsa = GetWedgeTexCoordStorageAttribute(m);
    for (auto& sf : shell.face) {
        CoordStorage target;
        // if the face is hole-filling the target geometry is its own shape, scaled according to the target geometry
        if (sf.IsHoleFilling()) {
            double scale = 1.0;
            if (targetGeometry == Model)
                scale = std::sqrt(avg3D / DistortionMetric::Area3D(sf));
            else if (targetGeometry == Texture)
                scale = std::sqrt(avgUV / DistortionMetric::Area3D(sf));
            else
                assert(0 && "Unexpected targetGeometry parameter value");
            target.P[0] = sf.P(0) * scale;
            target.P[1] = sf.P(1) * scale;
            target.P[2] = sf.P(2) * scale;
        } else {
            auto& mf = m.face[ia[sf]];
            if (targetGeometry == Model) {
                target.P[0] = mf.P(0);
                target.P[1] = mf.P(1);
                target.P[2] = mf.P(2);
            } else if (targetGeometry == Texture) {
                const Point2d& u0 = wtcsa[mf].tc[0].P();
                const Point2d& u1 = wtcsa[mf].tc[1].P();
                const Point2d& u2 = wtcsa[mf].tc[2].P();
                Point2d u10 = u1 - u0;
                Point2d u20 = u2 - u0;
                double area = std::abs(u10 ^ u20) / 2.0;
                if (area == 0) {
                    // zero uv area face, target the 3d shape scaled according to the target geometry
                    double scale = std::sqrt(avgUV / DistortionMetric::Area3D(mf));
                    target.P[0] = mf.P(0) * scale;
                    target.P[1] = mf.P(1) * scale;
                    target.P[2] = mf.P(2) * scale;
                } else {
                    // target shape is the original uv face
                    target.P[0] = Point3d{u0[0], u0[1], 0};
                    target.P[1] = Point3d{u1[0], u1[1], 0};
                    target.P[2] = Point3d{u2[0], u2[1], 0};
                }
            }
        }
        tsa[sf] = target;
    }
    return true;
}

void ChartOutlinesUV(Mesh& m, FaceGroup& chart, std::vector<std::vector<Point2f>> &outline2Vec)
{
    std::vector<std::vector<Point2d>> outVecDouble;
    ChartOutlinesUV(m, chart, outVecDouble);

    outline2Vec.clear();
    for (auto& vec : outVecDouble) {
        std::vector<Point2f> outline;
        for (auto& p : vec)
            outline.push_back(vcg::Point2f(p[0], p[1]));
        outline2Vec.push_back(outline);
    }
}

void ChartOutlinesUV(Mesh& m, FaceGroup& chart, std::vector<std::vector<Point2d>> &outline2Vec)
{
    assert(chart.numMerges == 0);

    struct FaceAdjacency {
        bool changed[3] = {false, false, false};
        Mesh::FacePointer FFp[3];
        char FFi[3];
    };

    /*
    auto CCIDh = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(m, CCIDh));

    std::unordered_map<Mesh::FacePointer,FaceAdjacency> savedAdj;
    for (auto fptr : chart.fpVec) {
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
    */

    std::unordered_set<Mesh::FacePointer> chartFaces;
    for (auto fptr : chart.fpVec)
        chartFaces.insert(fptr);

    std::unordered_map<Mesh::FacePointer,FaceAdjacency> savedAdj;
    for (auto fptr : chart.fpVec) {
        fptr->ClearV();
        for (int i = 0; i < 3; ++i) {
            if (chartFaces.count(fptr->FFp(i)) == 0) {
                savedAdj[fptr].changed[i] = true;
                savedAdj[fptr].FFp[i] = fptr->FFp(i);
                savedAdj[fptr].FFi[i] = fptr->FFi(i);
                fptr->FFp(i) = fptr;
                fptr->FFi(i) = i;
            }
        }
    }

    outline2Vec.clear();
    std::vector<Point2d> outline;

    for (auto fptr : chart.fpVec)
        fptr->ClearV();

    for (auto fptr : chart.fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (!fptr->IsV() && face::IsBorder(*fptr, i)) {
                face::Pos<Mesh::FaceType> p(fptr, i);
                face::Pos<Mesh::FaceType> startPos = p;
                assert(p.IsBorder());
                do {
                    assert(p.IsManifold());
                    p.F()->SetV();
                    //outline.push_back(Point2<ScalarType>(p.V()->P()));
                    outline.push_back(p.F()->WT(p.VInd()).P());
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

    std::size_t maxsz = 0;
    for (std::size_t i = 0; i < outline2Vec.size(); ++i) {
        maxsz = std::max(maxsz, outline2Vec[i].size());
    }
    if (maxsz == 0) {// fallback to uv box as outline
        vcg::Box2d box = chart.UVBox();
        outline.push_back(Point2d(box.min.X(), box.min.Y()));
        outline.push_back(Point2d(box.max.X(), box.min.Y()));
        outline.push_back(Point2d(box.max.X(), box.max.Y()));
        outline.push_back(Point2d(box.min.X(), box.max.Y()));
        outline2Vec.push_back(outline);
        outline.clear();
    }
}

bool CleanSmallComponents(Mesh& m, GraphHandle graph, TextureObjectHandle texObj, double areaThreshold)
{
    (void) texObj;

    auto CCIDh = GetConnectedComponentIDAttribute(m);
    tri::UpdateTopology<Mesh>::FaceFace(m);
    tri::UpdateTopology<Mesh>::VertexFace(m);
    tri::UpdateSelection<Mesh>::FaceClear(m);
    tri::UpdateFlags<Mesh>::FaceClearV(m);
    bool selected = false;
    for (auto& p : graph->charts) {
        ChartHandle chart = p.second;
        if (chart->AreaUV() < areaThreshold) {
            std::cout << "Chart " << chart->id << std::endl;
            std::unordered_set<RegionID> adj;
            for (auto ch : chart->adj)
                adj.insert(ch->id);
            std::vector<Mesh::FacePointer> visitVec;
            std::unordered_set<RegionID> visitedComponents;
            for (auto fptr : chart->fpVec) {
                for (int i = 0; i < 3; ++i) {
                    vcg::face::VFIterator<Mesh::FaceType> vfi(fptr->V(i));
                    while (!vfi.End()) {
                        auto fp = vfi.F();
                        if (!fp->IsV()) {
                            fp->SetV();
                            visitVec.push_back(fp);
                            visitedComponents.insert(CCIDh[fp]);
                        }
                        ++vfi;
                    }
                }
            }
            // reset visit flags
            for (auto fptr : visitVec)
                fptr->ClearV();
            // Only select if the visited components are 2, otherwise maintaining
            // reasonable uv coordinates is troublesome
            if (visitedComponents.size() <= 2) {
                // Check if the visited faces result in a surface of genus > 0
                Mesh probe;
                MeshFromFacePointers(visitVec, probe);
                tri::UpdateTopology<Mesh>::FaceFace(probe);
                int nonManifoldPrimitives = tri::Clean<Mesh>::CountNonManifoldEdgeFF(probe) + tri::Clean<Mesh>::CountNonManifoldVertexFF(probe);
                if (nonManifoldPrimitives || tri::Clean<Mesh>::MeshGenus(probe) > 0) {
                    // select the region to clean it later
                    for (auto& fp : visitVec)
                        fp->SetS();
                    selected = true;
                }
            }
        }
    }

    if (!selected)
        return false;

    // store in a vertex attribute the uv coordinates at the boundary
    auto uvattr = tri::Allocator<Mesh>::GetPerVertexAttribute<Point2d>(m, "uvattr");
    for (auto& v : m.vert)
        uvattr[v] = Point2d::Zero();
    for (auto& f : m.face) {
        if (f.IsS()) {
            for (int i = 0; i < 3; ++i) {
                if (!f.FFp(i)->IsS() && CCIDh[f] != CCIDh[f.FFp(i)]) {
                    // the adjacent face is not selected and the edge is a seam edge
                    uvattr[f.V0(i)] = f.cWT(i).P();
                    uvattr[f.V1(i)] = f.cWT((i+1)%3).P();
                }
            }
        }
    }

    for (auto& f : m.face)
        if (!f.IsD() && f.IsS())
            tri::Allocator<Mesh>::DeleteFace(m, f);

    tri::UpdateTopology<Mesh>::VertexFace(m);
    tri::Allocator<Mesh>::CompactEveryVector(m);
    tri::UpdateTopology<Mesh>::FaceFace(m);

    std::size_t oldFaces = m.face.size();

    tri::Hole<Mesh>::EarCuttingFill<tri::MinimumWeightEar<Mesh>>(m, 30, false);
    tri::Clean<Mesh>::RemoveUnreferencedVertex(m);
    tri::Clean<Mesh>::RemoveDuplicateVertex(m);

    tri::UpdateTopology<Mesh>::VertexFace(m);
    tri::Allocator<Mesh>::CompactEveryVector(m);
    tri::UpdateTopology<Mesh>::FaceFace(m);

    assert(m.face.size() > oldFaces);
    std::cout << "Added " << m.face.size() - oldFaces << " faces" << std::endl;

    for (std::size_t h = oldFaces; h < m.face.size(); ++h) {
        auto& f = m.face[h];
        assert(!f.IsD());
        for (int i = 0; i < 3; ++i) {
            if (uvattr[f.V(i)] != Point2d::Zero()) {
                f.WT(i).P() = uvattr[f.V(i)];
            } else {
                auto& fj = *(f.FFp(i));
                if (tri::Index(m, fj) < oldFaces) {
                    int j = f.FFi(i);
                    vcg::Point2d wuv = fj.WT((j+1)%3).P();
                    face::Pos<MeshFace> startPos(&f, i);
                    face::Pos<MeshFace> p = startPos;
                    do {
                        if (p.F()->WT(p.VInd()).P() == vcg::Point2d::Zero()) {
                            p.F()->WT(p.VInd()).P() = wuv;
                        }
                        p.NextE();
                    } while (p != startPos);
                }
            }
        }
    }

    tri::Allocator<Mesh>::DeletePerVertexAttribute(m, "uvattr");
    return true;
}



// FaceGroup class implementation
// ==============================

FaceGroup::FaceGroup(Mesh& m, const RegionID id_)
    : mesh{m},
      id{id_},
      fpVec{},
      adj{},
      numMerges{0},
      minMappedFaceValue{-1},
      maxMappedFaceValue{-1},
      borderUV{0},
      borderChanged{false}
{
}

void FaceGroup::Clear()
{
    id = INVALID_ID;
    fpVec.clear();
    adj.clear();
    numMerges = 0;
    minMappedFaceValue = -1;
    maxMappedFaceValue = -1;
    borderUV = 0.0;
    borderChanged = false;
}

void FaceGroup::AddFace(const Mesh::FacePointer fptr)
{
    fpVec.push_back(fptr);
    ParameterizationChanged();
}

double FaceGroup::OriginalAreaUV() const
{
    assert(HasWedgeTexCoordStorageAttribute(mesh));
    auto wtcsattr = GetWedgeTexCoordStorageAttribute(mesh);

    double areaUV = 0;
    for (auto fptr : fpVec) {
        const TexCoordStorage& tcs = wtcsattr[fptr];
        areaUV += std::abs((tcs.tc[1].P() - tcs.tc[0].P()) ^ (tcs.tc[2].P() - tcs.tc[0].P()));
    }
    return areaUV;
}

double FaceGroup::AreaUV() const
{
    double areaUV = 0;
    for (auto fptr : fpVec) areaUV += std::abs(DistortionMetric::AreaUV(*fptr));
    return areaUV;
}

double FaceGroup::Area3D() const
{
    double area3D = 0;
    for (auto fptr : fpVec) area3D += DistortionMetric::Area3D(*fptr);
    return area3D;
}

double FaceGroup::BorderUV() const
{
    if (borderChanged) UpdateBorder();
    return borderUV;
}

vcg::Box2d FaceGroup::UVBox() const
{
    vcg::Box2d box;
    for (auto fptr : fpVec) {
        box.Add(fptr->WT(0).P());
        box.Add(fptr->WT(1).P());
        box.Add(fptr->WT(2).P());
    }
    return box;
}

void FaceGroup::ParameterizationChanged()
{
    borderChanged = true;
}

Mesh::FacePointer FaceGroup::Fp()
{
    assert(!fpVec.empty()); return fpVec[0];
}

std::size_t FaceGroup::FN() const
{
    return fpVec.size();
}

std::size_t FaceGroup::NumAdj() const
{
    return adj.size();
}

void FaceGroup::MapDistortion(DistortionMetric::Type type, ParameterizationGeometry geometry, double areaScale)
{
    minMappedFaceValue = std::numeric_limits<float>::max();
    maxMappedFaceValue = std::numeric_limits<float>::lowest();
    for (auto fptr : fpVec) {
        if (type == DistortionMetric::Area) {
            fptr->Q() = DistortionMetric::AreaDistortion(mesh, *fptr, areaScale, geometry);
        }
        else if (type == DistortionMetric::Angle) {
            fptr->Q() = DistortionMetric::AngleDistortion(mesh, *fptr, geometry);
        }
        else assert(0 && "FaceGroup::MapDistortion");
        minMappedFaceValue = std::min(minMappedFaceValue, fptr->Q());
        maxMappedFaceValue = std::max(maxMappedFaceValue, fptr->Q());
    }
}

void FaceGroup::UpdateBorder() const
{
    auto CCIDh = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(mesh, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(mesh, CCIDh));

    borderUV = 0.0;
    for (auto fptr : fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (face::IsBorder(*fptr, i) || CCIDh[fptr] != CCIDh[fptr->FFp(i)]) {
                borderUV += EdgeLengthUV(*fptr, i);
            }
        }
    }
    borderChanged = false;
}


// MeshGraph class implementation
// ==============================

MeshGraph::MeshGraph(Mesh& m)
    : mesh{m}
{
}

void MeshGraph::MapDistortion(DistortionMetric::Type type, ParameterizationGeometry geometry)
{
    double areaScale;
    DistortionMetric::ComputeAreaScale(mesh, areaScale, geometry);
    /// FIXME absolute distortion should be an option parameter
    if (geometry == ParameterizationGeometry::Texture) areaScale = 1;
    for (auto& c : charts) c.second->MapDistortion(type, geometry, areaScale);

    // Map distortion to color
    std::pair<float, float> range = DistortionRange();
    for (auto& c : charts) {
        for (auto fptr : c.second->fpVec) {
            float q = fptr->Q();
            //float q = fptr->Q() / std::max(std::abs(range.first), std::abs(range.second));
            //float q = fptr->Q() / std::abs(range.first);
            if (q < 0) {
                //float v = 1.0f - (q / range.first);
                float v = 1.0f + (q / (-range.first));
                fptr->C().Import(Color4f{1.0f, v*v, v*v, 1.0f});
            } else {
                //float v = 1.0f - (q / range.second);
                float v = 1.0f - (q / range.second);
                fptr->C().Import(Color4f{v*v, v*v, 1.0f, 1.0f});
            }
        }
    }
}

std::pair<float,float> MeshGraph::DistortionRange() const
{
    std::pair<float,float> range = std::make_pair(std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest());
    for (const auto& c : charts) {
        range.first = std::min(c.second->minMappedFaceValue, range.first);
        range.second = std::max(c.second->maxMappedFaceValue, range.second);
    }
    return range;
}

std::shared_ptr<FaceGroup> MeshGraph::GetChart(RegionID i)
{
    //assert(charts.find(i) != charts.end() && "Chart does not exist");
    auto e = charts.find(i);
    if (e != charts.end()) return e->second;
    else return nullptr;
}

std::shared_ptr<FaceGroup> MeshGraph::GetChart_Insert(RegionID i)
{
    if (charts.find(i) == charts.end()) charts.insert(std::make_pair(i, std::make_shared<FaceGroup>(mesh, i)));
    return charts[i];
}

std::size_t MeshGraph::Count() const
{
    return charts.size();
}

int MeshGraph::MergeCount() const
{
    int n = 0;
    for (const auto& c : charts) n += c.second->numMerges;
    return n;

}

double MeshGraph::Area3D() const
{
    double area3D = 0.0f;
    for (const auto& c : charts) area3D += c.second->Area3D();
    return area3D;
}

double MeshGraph::AreaUV() const
{
    double areaUV = 0.0f;
    for (const auto& c : charts) areaUV += c.second->AreaUV();
    return areaUV;
}

double MeshGraph::BorderUV(float *meshBorderLengthUV, float *seamLengthUV)
{
    double borderUV = 0.0;
    for (const auto& c : charts) borderUV += c.second->BorderUV();
    return borderUV;
}

// GraphManager class implementation
// =================================

int GraphManager::Collapse_OK = 0;
int GraphManager::Collapse_ERR_DISCONNECTED = 1;
int GraphManager::Collapse_ERR_UNFEASIBLE = 2;

GraphManager::GraphManager(std::shared_ptr<MeshGraph> gptr, std::unique_ptr<EdgeWeightFunction> &&wfct)
    : g{gptr},
      wfct{std::move(wfct)}
{
    tri::UpdateTopology<Mesh>::FaceFace(g->mesh);
    for (auto elem : g->charts) {
        ChartHandle c1 = elem.second;
        for (ChartHandle c2 : c1->adj) {
            AddEdge(c1, c2);
        }
    }
}


std::shared_ptr<MeshGraph> GraphManager::Graph()
{
    return g;
}

bool GraphManager::Valid(std::pair<GraphManager::Edge,double> weightedEdge)
{
    auto e = edges.find(weightedEdge.first);
    return e != edges.end() && e->second == weightedEdge.second;
}

const std::pair<GraphManager::Edge,double>& GraphManager::PeekNextEdge()
{
    //assert(HasNextEdge());
    return queue.top();
}

void GraphManager::RemoveNextEdge()
{
    //assert(HasNextEdge());
    queue.pop();
}

bool GraphManager::HasNextEdge()
{
    constexpr double infinity = std::numeric_limits<double>::infinity();
    while (queue.size() > 0)
    {
        auto we = queue.top();
        if (Valid(we)) {
            if (we.second == infinity) { // The remaining edges cannot be collapsed due to infeasibility
                return false;
            }
            else {
                Mesh test;
                GraphManager::Edge e = we.first;
                BuildMeshFromFacePointers(test, {&(e.a->fpVec), &(e.b->fpVec)});
                if (Parameterizable(test)) {
                    return true;
                }
                else { // Cannot collapse the edge, reinsert with infinite weight
                    queue.pop();
                    we.second = infinity;
                    edges[we.first] = infinity;
                    queue.push(we);
                }
            }
        }
        else { // Invalid edge, throw it away
            queue.pop();
        }
    }
    return false;
}

/// TODO this is not robust if a client tries to collapse an arbitrary edge rather than the one returned
/// by *NextEdge() since feasibility is lazily evaluated (unfeasible edges may be present in the edge set)
ChartHandle GraphManager::Collapse(const GraphManager::Edge& e)
{
    ChartHandle c1, c2;
    if (e.a->FN() > e.b->FN()) {
        c1 = e.a; c2 = e.b;
    }
    else {
        c1 = e.b; c2 = e.a;
    }

    Merge(c1, c2);
    return c1;
}

bool GraphManager::AddEdge(ChartHandle c1, ChartHandle c2, bool replace)
{
     GraphManager::Edge e{c1, c2};

     if (replace) edges.erase(e);

     if (edges.find(e) == edges.end()) {
         auto weightedEdge = std::make_pair(e, (*wfct)(e));
         edges.insert(weightedEdge);
         queue.push(weightedEdge);
         return true;
     }
     else {
         return false;
     }
}

void GraphManager::Merge(ChartHandle c1, ChartHandle c2)
{
    assert(edges.find(GraphManager::Edge{c1,c2}) != edges.end());

    auto CCIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(g->mesh, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(g->mesh, CCIDh));

    // Update ID and add faces
    for (auto fp : c2->fpVec) {
        CCIDh[fp] = c1->id;
        c1->AddFace(fp);
    }

    // Update adjacencies - Note that obsolete edges remain in the priority queue
    c1->adj.erase(c2);
    for (auto cn : c2->adj) {
        GraphManager::Edge e2{c2, cn};
        edges.erase(e2);
        if (cn != c1) { // c1 is now (if it wasn't already) adjacent to cn
            cn->adj.erase(c2);
            cn->adj.insert(c1);
            c1->adj.insert(cn);

            // Manage queue and edge map state
            AddEdge(c1, cn, true); // replace edge
            /*
            GraphManager::Edge e1{c1, cn};
            edges.erase(e1); // delete the edge if it already exists so it can be reweighted
            auto weighted = std::make_pair(e1, wfct(e1));
            edges.insert(weighted);
            queue.push(weighted);
            */
        }
    }

    c1->numMerges += c2->numMerges + 1;

    g->charts.erase(c2->id);
}

void GraphManager::Split(const RegionID id, std::vector<ChartHandle>& splitCharts)
{
    assert(g->charts.count(id) == 1);

    assert(HasConnectedComponentIDAttribute(g->mesh));
    auto CCIDh = GetConnectedComponentIDAttribute(g->mesh);

    assert(HasInitialConnectedComponentIDAttribute(g->mesh));
    auto ICCh  = GetInitialConnectedComponentIDAttribute(g->mesh);

    ChartHandle chart = g->GetChart(id);
    if (chart->numMerges == 0) {
        splitCharts.push_back(chart);
        return;
    }

    // detach chart from the graph
    g->charts.erase(chart->id);
    for (auto cn : chart->adj) {
        cn->adj.erase(chart); // remove reference in the neighbor
        edges.erase(GraphManager::Edge{chart, cn}); // remove edge from support list
    }

    // restore original ids
    for (auto fptr : chart->fpVec) {
        CCIDh[fptr] = ICCh[fptr];
    }

    std::set<ChartHandle> newCharts;
    // readd faces to the graph, creating the required nodes
    for (auto fptr : chart->fpVec) {
        RegionID regionId = CCIDh[fptr];
        auto ch = g->GetChart_Insert(regionId);
        newCharts.insert(ch);
        ch->AddFace(fptr);
        for (int i = 0; i < fptr->VN(); ++i) {
            RegionID adjId = CCIDh[fptr->FFp(i)];
            if (regionId != adjId) {
                // force symmetry
                (ch->adj).insert(g->GetChart_Insert(adjId));
                (g->GetChart_Insert(adjId)->adj).insert(ch);
            }
        }
    }

    // reset merge count and insert new edges in the support list
    for (auto c1 : newCharts) {
        assert(c1->fpVec.size() > 0);
        c1->numMerges = 0;
        c1->id = CCIDh[c1->fpVec[0]];
        splitCharts.push_back(c1);
        for (auto c2 : c1->adj) {
            AddEdge(c1, c2, true);
        }
    }
}

int GraphManager::CloseMacroRegions(std::size_t minRegionSize)
{
    int mergeCount = 0;

    auto& regions = g->charts;

    std::unordered_map<RegionID, std::vector<RegionID>> mergeLists;
    std::unordered_map<RegionID, RegionID> invertedIndex;

    for (const auto& entry : regions) {
        auto chart = entry.second;
        if (invertedIndex.count(chart->id) == 1) continue; // skip if this region is already going to be merged to something else
        for (auto& adjRegion : chart->adj)
            if (adjRegion->NumAdj() == 1 && adjRegion->FN() < minRegionSize) {
                assert(invertedIndex.count(adjRegion->id) == 0);
                mergeLists[chart->id].push_back(adjRegion->id);
                invertedIndex[adjRegion->id] = chart->id;
            }
    }

    for (auto& entry : mergeLists) {
        Mesh probe;
        std::vector<std::vector<Mesh::FacePointer>* > fpVecp;
        fpVecp.reserve(entry.second.size()+1);
        fpVecp.push_back(&(regions[entry.first]->fpVec));
        for (auto id : entry.second) {
            fpVecp.push_back(&(regions[id]->fpVec));
        }
        BuildMeshFromFacePointers(probe, fpVecp);
        if (Parameterizable(probe)) {
            std::cout << "Merging " << entry.second.size() << " islands from macro region " << entry.first << std::endl;
            for (auto id : entry.second) {
                // don't bother checking feasibility of each merge since we already know that the union of all the charts is
                Merge(regions[entry.first], regions[id]);
                mergeCount++;
            }
        } else {
            //tri::io::ExporterOBJ<PMesh>::Save(probe, "fail.obj", tri::io::Mask::IOM_WEDGTEXCOORD);
        }
    }

    return mergeCount;
}
