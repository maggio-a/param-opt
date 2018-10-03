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
#include "utils.h"

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


inline Polyline2 BuildPolyline2(const std::vector<std::size_t> &vfi, const std::vector<int> &vvi, const Mesh& shell2D)
{
    Polyline2 polyline;
    polyline.reserve(vfi.size());
    for (std::size_t i = 0; i < vfi.size(); ++i) {
        const vcg::Point3d& p = shell2D.face[vfi[i]].cP(vvi[i]);
        ensure_condition(p.Z() == 0.0);
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
                ensure_condition(p.IsBorder());
                do {
                    ensure_condition(p.IsManifold());
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

#include <vcg/complex/algorithms/isotropic_remeshing.h>

void CloseMeshHoles(Mesh& shell)
{
    int startFN = shell.FN();

    // Get border info
    ensure_condition(HasBoundaryInfoAttribute(shell));
    BoundaryInfo& info = GetBoundaryInfoAttribute(shell)();

    //vcg::tri::io::ExporterPLY<Mesh>::Save(shell, "original.ply", tri::io::Mask::IOM_FACECOLOR);
    // Leave only the longest boundary
    tri::UpdateFlags<Mesh>::FaceClearS(shell);
    ensure_condition(info.vBoundaryFaces.size() > 0 && "Mesh has no boundaries");
    if (info.vBoundaryFaces.size() > 1) {
        std::size_t k = info.LongestBoundary();
        // select all the boundary faces
        for (std::size_t i = 0; i < info.vBoundaryFaces.size(); ++i) {
            if (i == k) continue;
            for (auto j : info.vBoundaryFaces[i]) {
                ensure_condition(face::IsBorder(shell.face[j], 0) || face::IsBorder(shell.face[j], 1) || face::IsBorder(shell.face[j], 2));
                shell.face[j].SetS();
            }
        }
        tri::Hole<Mesh>::EarCuttingFill<tri::MinimumWeightEar<Mesh>>(shell, shell.FN(), true);
    }

    tri::Clean<Mesh>::RemoveZeroAreaFace(shell);
    tri::Allocator<Mesh>::CompactEveryVector(shell);

    // Compute border info
    double length = 0;
    int count = 0;
    for (std::size_t i = 0; i < info.vBoundaryFaces.size(); ++i) {
        if (i == info.LongestBoundary())
            continue;
        length += info.vBoundaryLength[i];
        count += (int) info.vBoundarySize[i];
    }

    // select the hole-filling faces that need to be remeshed
    tri::UpdateFlags<Mesh>::FaceClearS(shell);
    for (auto& sf : shell.face) {
        if (sf.IsHoleFilling()) sf.SetS();
    }

    /*
    int startFN = 0;
    int max_i  = 0;
    for (auto& sf : shell.face) {
        if (sf.IsMesh()) {
            assert(ia[sf] != -1);
            startFN++;
            max_i = std::max(max_i, int(tri::Index<Mesh>(shell, sf)));
        }
    }
    */

    // Remesh filled hole
    //ColorFace(shell);
    //vcg::tri::io::ExporterPLY<Mesh>::Save(shell, "original_closed.ply", tri::io::Mask::IOM_FACECOLOR);
    IsotropicRemeshing<Mesh>::Params params;
    //params.SetTargetLen(2.0*(totalBorderLen / totalBorderFaces));
    //params.SetTargetLen(length / count);
    params.SetTargetLen(length / (2 * count));
    params.SetFeatureAngleDeg(30);
    params.selectedOnly = true;
    //params.splitFlag = false;
    //params.swapFlag = false;
    //params.collapseFlag = false;
    params.smoothFlag = false;
    params.projectFlag = false;
    params.iter = 3;
    int iter = 0;
    do {
        IsotropicRemeshing<Mesh>::Do(shell, params);
        iter += params.iter;
    } while (params.stat.collapseNum + params.stat.flipNum + params.stat.splitNum > 0 && iter < 3);

    tri::Allocator<Mesh>::CompactEveryVector(shell);

    auto ia = GetFaceIndexAttribute(shell);
    tri::UpdateFlags<Mesh>::FaceClearS(shell);
    for (auto& f: shell.face) {
        if (int(tri::Index(shell, f)) >= startFN) {
            f.SetHoleFilling();
            ia[f] = -1;
        }
    }

    tri::UpdateTopology<Mesh>::FaceFace(shell);
    tri::UpdateTopology<Mesh>::VertexFace(shell);
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

// TODO this method should CALL ClearHoleFillingFaces(shell, true, true) otherwise it's unsafe
void SyncShellWithModel(Mesh& shell, Mesh& baseMesh)
{
    auto ia = GetFaceIndexAttribute(shell);
    for (auto& sf : shell.face) {
        ensure_condition(sf.IsMesh());
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

        ensure_condition(points.size() >= szoutline);

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
                ensure_condition(0 && "Unexpected targetGeometry parameter value");
            ensure_condition(scale > 0);
            target.P[0] = fi->P(0) * scale;
            target.P[1] = fi->P(1) * scale;
            target.P[2] = fi->P(2) * scale;
            tsa[fi] = target;

            Stabilize(tsa[fi]);

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
    ensure_condition(points.size() >= szoutline);
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

        Stabilize(tsa[fi]);

        fi++;
    }

    tri::UpdateTopology<Mesh>::FaceFace(shell);
    tri::UpdateTopology<Mesh>::VertexFace(shell);
}

#include <vcg/complex/algorithms/cut_tree.h>
#include <vcg/complex/algorithms/curve_on_manifold.h>
#include <vcg/complex/algorithms/crease_cut.h>

bool BuildShell(Mesh& shell, FaceGroup& fg, ParameterizationGeometry targetGeometry, bool useExistingUV)
{
    Mesh& m = fg.mesh;

    CopyToMesh(fg, shell);

    tri::Clean<Mesh>::RemoveDuplicateVertex(shell);

    tri::UpdateBounding<Mesh>::Box(shell);

    tri::UpdateTopology<Mesh>::FaceFace(shell);

    int splitCount = tri::Clean<Mesh>::SplitNonManifoldVertex(shell, 0.15);
    if (splitCount > 0) {
        std::cout << "Mesh was not vertex-manifold, " << splitCount << " vertices split" << std::endl;
    }
    auto ia = GetFaceIndexAttribute(shell);

    bool init = true;

    if (useExistingUV) {
        // copy wedge texture data from the input mesh
        for (auto& sf : shell.face) {
            auto& f = m.face[ia[sf]];
            for (int i = 0; i < 3; ++i) {
                sf.WT(i).P() = f.WT(i).P();
            }
        }
        // split vertices at seams
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
        if (!tri::Clean<Mesh>::IsCoherentlyOrientedMesh(shell)) {
            std::cout << "Attempting to reorient faces" << std::endl;
            bool p1, p2;
            tri::Clean<Mesh>::OrientCoherentlyMesh(shell, p1, p2);
            if (!p1)
                std::cout << "Mesh is not coherently oriented" << std::endl;
            if (!p2) {
                std::cout << "Mesh is non-orientable" << std::endl;
                tri::io::Exporter<Mesh>::Save(shell, "non_orientable.obj", tri::io::Mask::IOM_ALL);
            }
        }

        if (!Parameterizable(shell)) {
            Mesh cut;
            std::srand(0);
            tri::CutTree<Mesh> ct(shell);
            ct.Build(cut, rand() % shell.FN());
            tri::CoM<Mesh> cc(shell);
            cc.Init();
            bool ok = cc.TagFaceEdgeSelWithPolyLine(cut);
            ensure_condition(ok && "topological cut failed");
            //tri::io::Exporter<Mesh>::Save(shell, "to_cut.obj", tri::io::Mask::IOM_ALL);
            tri::CutMeshAlongSelectedFaceEdges(shell);
            tri::UpdateTopology<Mesh>::FaceFace(shell);
        }

        // Compute the initial configuration (Tutte's parameterization)
        ComputeBoundaryInfo(shell);
        CloseMeshHoles(shell);

        bool solved = false;

        if (!solved) {
            UniformSolver<Mesh> solver(shell);
            solved = solver.Solve();
            if (!solved)
                tri::io::Exporter<Mesh>::Save(shell, "error_tutte.obj", tri::io::Mask::IOM_ALL);
        }

        if (!solved) {
            UniformSolver<Mesh> solver(shell);
            solver.SetBoundaryMapProportional();
            solved = solver.Solve();
            if (!solved)
                tri::io::Exporter<Mesh>::Save(shell, "error_tutte_prop.obj", tri::io::Mask::IOM_ALL);
        }

        if (!solved) {
            MeanValueSolver<Mesh> mvs(shell);
            solved = mvs.Solve();
            if (!solved)
                tri::io::Exporter<Mesh>::Save(shell, "error_mvs.obj", tri::io::Mask::IOM_ALL);
        }

        if (!solved) {
            MeanValueSolver<Mesh> mvs(shell);
            mvs.UseCotangentWeights();
            solved = mvs.Solve();
            if (!solved)
                tri::io::Exporter<Mesh>::Save(shell, "error_cotan.obj", tri::io::Mask::IOM_ALL);
        }

        init = solved;
    }

    // if we reach this point, the shell is synced with its uv coordinates

    // Compute average areas
    auto psi = GetParameterizationScaleInfoAttribute(m);
    double avg3D = psi().surfaceArea / psi().numNonZero;
    double avgUV = psi().parameterArea / psi().numNonZero;

    double targetArea = 0;

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
                ensure_condition(0 && "Unexpected targetGeometry parameter value");
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
                //Point2d u10 = u1 - u0;
                //Point2d u20 = u2 - u0;
                Point2d u10;
                Point2d u20;
                LocalIsometry(u1 - u0, u2 - u0, u10, u20);
                double area = std::abs(u10 ^ u20) / 2.0;

                if (area == 0) {
                    // zero uv area face, target the 3d shape scaled according to the target geometry
                    double scale = std::sqrt(avgUV / DistortionMetric::Area3D(mf));
                    target.P[0] = mf.P(0) * scale;
                    target.P[1] = mf.P(1) * scale;
                    target.P[2] = mf.P(2) * scale;
                } else {
                    // Compute the isometry to the model coordinates and scale them to match the texture area
                    Point2d x10;
                    Point2d x20;
                    LocalIsometry(mf.P(1) - mf.P(0), mf.P(2) - mf.P(0), x10, x20);

                    //double areaModel = std::abs(x10 ^ x20) / 2.0;
                    //x10 *= std::sqrt(area / areaModel);
                    //x20 *= std::sqrt(area / areaModel);

                    // scale using the average scaling factor rather than the one of
                    // the triangle, since heavily distorted triangle are likely to
                    // have areas that are too small
                    x10 *= psi().scale;
                    x20 *= psi().scale;

                    // compute the singular values of the transformation matrix, s2 > s1
                    // ref for the formula: smith&schaefer 2015 bijective,
                    Eigen::Matrix2d phi = ComputeTransformationMatrix(x10, x20, u10, u20);
                    double bcplus  = std::pow(phi(0, 1) + phi(1, 0), 2.0);
                    double bcminus = std::pow(phi(0, 1) - phi(1, 0), 2.0);
                    double adplus  = std::pow(phi(0, 0) + phi(1, 1), 2.0);
                    double adminus = std::pow(phi(0, 0) - phi(1, 1), 2.0);
                    double s_min = 0.5 * std::abs(std::sqrt(bcplus + adminus) - std::sqrt(bcminus + adplus));
                    double s_max = 0.5 * (std::sqrt(bcplus + adminus) + std::sqrt(bcminus + adplus));

                    double interpolationFactor = 1.0 - (s_min / s_max);
                    ensure_condition(interpolationFactor > 0);
                    ensure_condition(interpolationFactor <= 1);
                    ensure_condition(std::isfinite(interpolationFactor));

                    Point2d v10 = (1.0 - interpolationFactor) * u10 + (interpolationFactor) * x10;
                    Point2d v20 = (1.0 - interpolationFactor) * u20 + (interpolationFactor) * x20;

                    target.P[0] = Point3d(0, 0, 0);
                    target.P[1] = Point3d(v10[0], v10[1], 0);
                    target.P[2] = Point3d(v20[0], v20[1], 0);
                }
            }
        }
        tsa[sf] = target;
        //Point2d u10 = Point2d(target.P[1].X() - target.P[0].X(), target.P[1].Y() - target.P[0].Y());
        //Point2d u20 = Point2d(target.P[2].X() - target.P[0].X(), target.P[2].Y() - target.P[0].Y());
        //targetArea += std::abs(u10 ^ u20) / 2.0;
        targetArea += ((target.P[1] - target.P[0]) ^ (target.P[2] - target.P[0])).Norm() / 2.0;
        /*
        std::cout << target.P[0].X() << "  " << target.P[0].Y() << "  " << target.P[0].Z() << std::endl;
        std::cout << target.P[1].X() << "  " << target.P[1].Y() << "  " << target.P[1].Z() << std::endl;
        std::cout << target.P[2].X() << "  " << target.P[2].Y() << "  " << target.P[2].Z() << std::endl;
        std::cout << "------------ " << (sf.IsHoleFilling() == true) << std::endl;
        */
        ensure_condition(std::isfinite(targetArea));
    }

    for (auto& sf : shell.face) {
        Stabilize(tsa[sf]);
    }

    double shellUvArea = 0;
    for (auto& f : shell.face) {
        double area = ((f.V(1)->T().P() - f.V(0)->T().P()) ^ (f.V(2)->T().P() - f.V(0)->T().P())) / 2.0;
        shellUvArea += std::abs(area);
    }

    // Scale the shell size to match the target shape area
    double areaScaleFactor = std::sqrt(targetArea / shellUvArea);
    std::cout << targetArea << " " << shellUvArea << std::endl;
    ensure_condition(areaScaleFactor > 0);
    for (auto& v : shell.vert)
        v.T().P() *= areaScaleFactor;
    SyncShellWithUV(shell);

    return init;
}

/* This function is used to 'correct' degenerate triangles (too small and/or
 * slivers), by assigning to the CoordStorage ref the shape of a equiareal
 * triangle of comparable area. This should prevent numerical issues during
 * the optimization process */
void Stabilize(CoordStorage& cs)
{
    /*
    Point3d p10 = cs.P[1] - cs.P[0];
    Point3d p20 = cs.P[2] - cs.P[0];
    double targetArea = (p10 ^ p20).Norm() / 2.0;
    double q = QualityRadii(Point3d::Zero(), Point3d::Zero() + p10, Point3d::Zero() + p20);
    if (q < 0.1) {
        p10 = Point3d(1, 0, 0);
        p20 = Point3d(std::cos(M_PI / 3.0), std::sin(M_PI / 3.0), 0);
        p10 *= std::sqrt(targetArea);
        p20 *= std::sqrt(targetArea);
    }
    cs.P[1] = cs.P[0] + p10;
    cs.P[2] = cs.P[0] + p20;
    */
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
    ensure_condition(chart.numMerges == 0);

    outline2Vec.clear();
    std::vector<Point2d> outline;

    for (auto fptr : chart.fpVec)
        fptr->ClearV();

    for (auto fptr : chart.fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (!fptr->IsV() && face::IsBorder(*fptr, i)) {
                face::Pos<Mesh::FaceType> p(fptr, i);
                face::Pos<Mesh::FaceType> startPos = p;
                ensure_condition(p.IsBorder());
                do {
                    ensure_condition(p.IsManifold());
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

    std::size_t maxsz = 0;
    for (std::size_t i = 0; i < outline2Vec.size(); ++i) {
        maxsz = std::max(maxsz, outline2Vec[i].size());
    }
    if (maxsz == 0) {// fallback to uv box as outline
        std::cout << "[LOG] Falling back to UVBox as outline for chart " << chart.id << std::endl;
        vcg::Box2d box = chart.UVBox();
        outline.push_back(Point2d(box.min.X(), box.min.Y()));
        outline.push_back(Point2d(box.max.X(), box.min.Y()));
        outline.push_back(Point2d(box.max.X(), box.max.Y()));
        outline.push_back(Point2d(box.min.X(), box.max.Y()));
        outline2Vec.push_back(outline);
        outline.clear();
    }
}

/*
 * This function tries to eliminate small handles that can be identified by regions
 * of genus non-zero surrounding very small charts.
 * The logic is
 *   1. Iterate over small charts (those with uv area below the areaThreshold parameter
 *     1.1 Visit the faces of the chart and those adjacent to any visited face/vertex
 *     1.2 If the visited region has genus>0 and 'the texture coordinates can be
 *         recovered' (*), select the region
 *   2. Delete all the selected faces
 *   3. Close the holes that are left (by filling ears)
 *   4. Restore the texture coordinates around the closed hole
 *
 * (*) Recovering texture coordinates
 * To ensure that valid texture coordinates are kept around the filled area, check
 * that at most one chart is not fully selected. This guarantees that no seams are
 * crossed by the selection.
 * Note that if the selection only touches a seam, we store the wedge texcoords
 * of the face to be deleted on the vertices (as attributes) to easily retrieve
 * them after the face is deleted.
 * */
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
            /*
            std::unordered_set<RegionID> adj;
            for (auto ch : chart->adj)
                adj.insert(ch->id);
                */
            std::vector<Mesh::FacePointer> visitVec;
            std::unordered_set<RegionID> visitedComponents;

            // I do the visit in two passes, first i select faces that are edge-adjacent to
            // the seed chart, then I expand the visit around vertices only to the region
            // already visited

            for (auto fptr : chart->fpVec) {
                fptr->SetV();
                visitVec.push_back(fptr);
                visitedComponents.insert(CCIDh[fptr]);
            }
            for (auto fptr : chart->fpVec) {
                for (int i = 0; i < 3; ++i) {
                    Mesh::FacePointer fp = fptr->FFp(i);
                    if (!fp->IsV()) {
                        fp->SetV();
                        ensure_condition(std::find(visitVec.begin(), visitVec.end(), fp) == visitVec.end());
                        visitVec.push_back(fp);
                        visitedComponents.insert(CCIDh[fp]);
                    }
                }
            }

            for (auto fptr : chart->fpVec) {
                for (int i = 0; i < 3; ++i) {
                    vcg::face::VFIterator<Mesh::FaceType> vfi(fptr->V(i));
                    while (!vfi.End()) {
                        auto fp = vfi.F();
                        if (!fp->IsV() && (visitedComponents.count(CCIDh[fp]) > 0)) {
                            fp->SetV();
                            ensure_condition(std::find(visitVec.begin(), visitVec.end(), fp) == visitVec.end());
                            visitVec.push_back(fp);
                            visitedComponents.insert(CCIDh[fp]);
                        }
                        ++vfi;
                    }
                }
            }

            int numNotFullyVisited = 0;
            for (auto id : visitedComponents) {
                for (auto fptr : graph->GetChart(id)->fpVec) {
                    if (!fptr->IsV()) {
                        numNotFullyVisited++;
                        break;
                    }
                }
            }

            // Check that only one region touches seams, otherwise the wedge
            // texcoords will not agree
            std::unordered_set<RegionID> boundaryId;
            for (auto fptr : visitVec) {
                for (int i = 0; i < 3; ++i) {
                    if (CCIDh[fptr] != CCIDh[fptr->FFp(i)] && !fptr->FFp(i)->IsV()) {
                        boundaryId.insert(CCIDh[fptr]);
                    }
                }
            }

            // Only select the region if at most one region was not fully visited
            // this ensures that either a (small) component disappears, or we are
            // able to recover meaningful texture coordinates at the boundary of
            // the selection.
            if (boundaryId.size() <= 1 && numNotFullyVisited <= 1) {
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

            // reset visit flags
            for (auto fptr : visitVec)
                fptr->ClearV();
        }
    }

    if (!selected)
        return false;


    for (auto& f : m.face)
        if (!f.IsD() && f.IsS())
           f.C() = vcg::Color4b::Cyan;
       else
           f.C() = vcg::Color4b::DarkGray;

    tri::io::Exporter<Mesh>::Save(m, "color.obj", io::Mask::IOM_ALL);

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

    ensure_condition(m.face.size() > oldFaces);
    std::cout << "Added " << m.face.size() - oldFaces << " faces" << std::endl;

    for (std::size_t h = oldFaces; h < m.face.size(); ++h) {
        auto& f = m.face[h];
        ensure_condition(!f.IsD());
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

    tri::io::Exporter<Mesh>::Save(m, "fixed.obj", io::Mask::IOM_ALL);

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
      cache{},
      minMappedFaceValue{-1},
      maxMappedFaceValue{-1}
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
    cache = {};
    dirty = false;
}

void FaceGroup::UpdateCache() const
{
    double areaUV = 0;
    double area3D = 0;
    vcg::Point3d weightedSumNormal = vcg::Point3d::Zero();
    for (auto fptr : fpVec) {
        areaUV += std::abs(DistortionMetric::AreaUV(*fptr));
        area3D += DistortionMetric::Area3D(*fptr);
        weightedSumNormal += (fptr->P(1) - fptr->P(0)) ^ (fptr->P(2) ^ fptr->P(0));
    }

    // TODO this does not take cuts into account...
    auto CCIDh = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(mesh, "ConnectedComponentID");
    ensure_condition(tri::Allocator<Mesh>::IsValidHandle<RegionID>(mesh, CCIDh));

    double borderUV = 0.0;
    for (auto fptr : fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (face::IsBorder(*fptr, i) || CCIDh[fptr] != CCIDh[fptr->FFp(i)]) {
                borderUV += EdgeLengthUV(*fptr, i);
            }
        }
    }

    double border3D = 0.0;
    for (auto fptr : fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (face::IsBorder(*fptr, i) || CCIDh[fptr] != CCIDh[fptr->FFp(i)]) {
                border3D += EdgeLength(*fptr, i);
            }
        }
    }

    cache.area3D = area3D;
    cache.areaUV = areaUV;
    cache.borderUV = borderUV;
    cache.border3D = border3D;
    cache.weightedSumNormal = weightedSumNormal;

    dirty = false;
}

vcg::Point3d FaceGroup::AverageNormal() const
{
    if (dirty)
        UpdateCache();
    vcg::Point3d avgN = ((cache.weightedSumNormal) / (2.0 * cache.area3D));
    return avgN.Normalize();
}

void FaceGroup::AddFace(const Mesh::FacePointer fptr)
{
    fpVec.push_back(fptr);
    dirty = true;
}

double FaceGroup::OriginalAreaUV() const
{
    ensure_condition(HasWedgeTexCoordStorageAttribute(mesh));
    auto wtcsattr = GetWedgeTexCoordStorageAttribute(mesh);

    double doubleAreaUV = 0;
    for (auto fptr : fpVec) {
        const TexCoordStorage& tcs = wtcsattr[fptr];
        doubleAreaUV += std::abs((tcs.tc[1].P() - tcs.tc[0].P()) ^ (tcs.tc[2].P() - tcs.tc[0].P()));
    }
    return 0.5 * doubleAreaUV;
}

double FaceGroup::AreaUV() const
{
    if (dirty)
        UpdateCache();
    return cache.areaUV;
}

double FaceGroup::Area3D() const
{
    if (dirty)
        UpdateCache();
    return cache.area3D;
}

double FaceGroup::BorderUV() const
{
    if (dirty)
        UpdateCache();
    return cache.borderUV;
}

double FaceGroup::Border3D() const
{
    if (dirty)
        UpdateCache();
    return cache.border3D;
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
    dirty = true;
}

Mesh::FacePointer FaceGroup::Fp()
{
    ensure_condition(!fpVec.empty()); return fpVec[0];
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
        else ensure_condition(0 && "FaceGroup::MapDistortion");
        minMappedFaceValue = std::min(minMappedFaceValue, fptr->Q());
        maxMappedFaceValue = std::max(maxMappedFaceValue, fptr->Q());
    }
}

void FaceGroup::UpdateBorder() const
{
    if (dirty)
        UpdateCache();
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
    //ensure_condition(charts.find(i) != charts.end() && "Chart does not exist");
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
    double area3D = 0;
    for (const auto& c : charts) area3D += c.second->Area3D();
    return area3D;
}

double MeshGraph::MappedFraction() const
{
    double area3D = 0;
    double mappedArea3D = 0;
    for (const auto& c : charts) {
        area3D += c.second->Area3D();
        if (c.second->AreaUV() > 0)
            mappedArea3D += c.second->Area3D();
    }
    return mappedArea3D / area3D;
}

double MeshGraph::AreaUV() const
{
    double areaUV = 0;
    for (const auto& c : charts) areaUV += c.second->AreaUV();
    return areaUV;
}

double MeshGraph::BorderUV() const
{
    double borderUV = 0;
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
    //ensure_condition(HasNextEdge());
    return queue.top();
}

void GraphManager::RemoveNextEdge()
{
    //ensure_condition(HasNextEdge());
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
                GraphManager::Edge e = we.first;
                std::vector<Mesh::FacePointer> fpv;
                fpv.reserve(e.a->FN() + e.b->FN());
                fpv.insert(fpv.end(), e.a->fpVec.begin(), e.a->fpVec.end());
                fpv.insert(fpv.end(), e.b->fpVec.begin(), e.b->fpVec.end());
                Mesh test;
                MeshFromFacePointers(fpv, test);
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

bool GraphManager::ExistsEdge(GraphManager::Edge e)
{
    return edges.count(e) > 0;
}

double GraphManager::EdgeWeight(GraphManager::Edge e)
{
    if (ExistsEdge(e))
        return (*wfct)(e);
    else
        return std::numeric_limits<double>::infinity();
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
    ensure_condition(edges.find(GraphManager::Edge{c1,c2}) != edges.end());

    auto CCIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(g->mesh, "ConnectedComponentID");
    ensure_condition(tri::Allocator<Mesh>::IsValidHandle<RegionID>(g->mesh, CCIDh));

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
    ensure_condition(g->charts.count(id) == 1);

    ensure_condition(HasConnectedComponentIDAttribute(g->mesh));
    auto CCIDh = GetConnectedComponentIDAttribute(g->mesh);

    ensure_condition(HasInitialConnectedComponentIDAttribute(g->mesh));
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
        ensure_condition(c1->fpVec.size() > 0);
        c1->numMerges = 0;
        c1->id = CCIDh[c1->fpVec[0]];
        splitCharts.push_back(c1);
        for (auto c2 : c1->adj) {
            AddEdge(c1, c2, true);
        }
    }
}

int GraphManager::CloseMacroRegions(double areaThreshold)
{
    ensure_condition(areaThreshold > 0 && areaThreshold < 1);

    int mergeCount = 0;

    auto& regions = g->charts;

    double thresholdValue = Graph()->Area3D() * areaThreshold;

    std::unordered_map<RegionID, std::vector<RegionID>> mergeLists;
    std::unordered_map<RegionID, RegionID> invertedIndex;

    for (const auto& entry : regions) {
        auto chart = entry.second;
        if (invertedIndex.count(chart->id) == 1) continue; // skip if this region is already going to be merged to something else
        for (auto& adjRegion : chart->adj) {
            if ((adjRegion->NumAdj() == 1) && (adjRegion->Area3D() < thresholdValue)) {
                ensure_condition(invertedIndex.count(adjRegion->id) == 0);
                mergeLists[chart->id].push_back(adjRegion->id);
                invertedIndex[adjRegion->id] = chart->id;
            }
        }
    }

    for (auto& entry : mergeLists) {

        std::vector<Mesh::FacePointer> fpv;

        fpv.insert(fpv.end(), regions[entry.first]->fpVec.begin(), regions[entry.first]->fpVec.end());
        for (auto id : entry.second)
            fpv.insert(fpv.end(), regions[id]->fpVec.begin(), regions[id]->fpVec.end());

        Mesh probe;
        MeshFromFacePointers(fpv, probe);

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
