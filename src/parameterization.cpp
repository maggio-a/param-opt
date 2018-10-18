#include "parameterization.h"
#include "mesh_graph.h"
#include "energy.h"
#include "iterative.h"
#include "mesh_utils.h"
#include "math_utils.h"
#include "timer.h"
#include "uniform_solver.h"
#include "mean_value_param.h"

#include "polygon2_triangulator.h"

#include <memory>

#include <wrap/io_trimesh/export.h>

#include <vcg/complex/algorithms/crease_cut.h>
#include <vcg/complex/algorithms/update/color.h>

#include <vcg/complex/algorithms/cut_tree.h>
#include <vcg/complex/algorithms/curve_on_manifold.h>

#include <Eigen/Core>

#include <wrap/io_trimesh/export.h>


/* Shell-related functions
 * ======================= */


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

bool BuildShell(Mesh& shell, FaceGroup& fg, ParameterizationGeometry targetGeometry, bool useExistingUV)
{
    Mesh& m = fg.mesh;

    CopyToMesh(fg, shell);

    tri::Clean<Mesh>::RemoveDuplicateVertex(shell);

    tri::UpdateBounding<Mesh>::Box(shell);

    tri::UpdateTopology<Mesh>::FaceFace(shell);

    int splitCount;
    while ((splitCount = tri::Clean<Mesh>::SplitNonManifoldVertex(shell, 0.15)) > 0) {
        std::cout << "Mesh was not vertex-manifold, " << splitCount << " vertices split" << std::endl;
        tri::Allocator<Mesh>::CompactEveryVector(shell);
    }
    //tri::io::Exporter<Mesh>::Save(shell, "shell.obj", tri::io::Mask::IOM_ALL);

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
            tri::CutMeshAlongSelectedFaceEdges(shell);
            tri::UpdateTopology<Mesh>::FaceFace(shell);

            if (!Parameterizable(shell)) {
                tri::UpdateFlags<Mesh>::FaceClearFaceEdgeS(shell);
                Mesh::FacePointer fp = &shell.face[0];
                fp->SetFaceEdgeS(0);
                fp->SetFaceEdgeS(1);
                fp->FFp(0)->SetFaceEdgeS(fp->FFi(0));
                fp->FFp(1)->SetFaceEdgeS(fp->FFi(1));
                tri::CutMeshAlongSelectedFaceEdges(shell);
                tri::UpdateTopology<Mesh>::FaceFace(shell);
                ensure_condition(Parameterizable(shell));
            }
        }

        tri::io::Exporter<Mesh>::Save(shell, "shell_initial.obj", tri::io::Mask::IOM_ALL);

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

    tri::io::Exporter<Mesh>::Save(shell, "shell_closed.obj", tri::io::Mask::IOM_ALL);

    // Compute the target shapes
    auto tsa = GetTargetShapeAttribute(shell);
    auto wtcsa = GetWedgeTexCoordStorageAttribute(m);
    for (auto& sf : shell.face) {
        CoordStorage target;
        // if the face is hole-filling the target triangle is its own shape, scaled according to the target geometry
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
                /*
                Point2d v10, v20;
                LocalIsometry(mf.P(1) - mf.P(0), mf.P(2) - mf.P(0), v10, v20);

                target.P[0] = Point3d(0, 0, 0);
                target.P[1] = Point3d(v10.X(), v10.Y(), 0);
                target.P[2] = Point3d(v20.X(), v20.Y(), 0);
                */

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
                    ensure_condition(interpolationFactor >= 0);
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
        targetArea += ((target.P[1] - target.P[0]) ^ (target.P[2] - target.P[0])).Norm() / 2.0;
        ensure_condition(std::isfinite(targetArea));
        ensure_condition(targetArea > 0);
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
    ensure_condition(areaScaleFactor > 0);

    for (auto& v : shell.vert)
        v.T().P() *= areaScaleFactor;

    SyncShellWithUV(shell);

    return init;
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





/* Constructor, destructor and utility methods
 * =========================================== */


ParameterizerObject::ParameterizerObject(ChartHandle c, ParameterizationStrategy strat)
    : shell{},
      baseMesh{c->mesh},
      chart{c},
      strategy(strat),
      energy{},
      opt{},
      needsRemeshing{false},
      iterationCount{0},
      gradientNormTolerance{1e-3},
      energyDiffTolerance{1e-6},
      state{ParameterizerState::OK}
{
    Reset();
}

ParameterizerObject::~ParameterizerObject()
{
    // empty
}

void ParameterizerObject::Reset()
{
    energy = nullptr;
    opt = nullptr;
    iterationCount = 0;
    shell.Clear();

    if (strategy.warmStart)
        std::cout << "ParameterizerObject: Using warm start" << std::endl;

    bool init = BuildShell(shell, *chart, strategy.geometry, strategy.warmStart);
    if (!init) {
        std::cout << "WARNING: Failed to initialize injective solution" << std::endl;
        state = ParameterizerState::InitializationFailed;
        return;
    }

    //ensure_condition(init && "Failed to initialize ParameterizerObject solution");

    if (strategy.scaffold)
        BuildScaffold(shell, strategy.geometry, baseMesh);

    if (strategy.padBoundaries == false)
        ClearHoleFillingFaces(shell, true, false);
    InitializeOptimizer();
    SyncShellWithUV(shell);
    MarkInitialSeamsAsFaux(shell, baseMesh);

    auto tsa = GetTargetShapeAttribute(shell);
    targetArea = 0;
    for (auto& sf : shell.face) {
        CoordStorage target = tsa[sf];
        // changing this because it seems more robust
        //Point2d u10 = Point2d(target.P[1].X() - target.P[0].X(), target.P[1].Y() - target.P[0].Y());
        //Point2d u20 = Point2d(target.P[2].X() - target.P[0].X(), target.P[2].Y() - target.P[0].Y());
        targetArea += ((target.P[1] - target.P[0]) ^ (target.P[2] - target.P[0])).Norm() / 2.0;
    }

//    std::cout << "ENERGY WHEN CONSTRUCTED == " << energy->E() << std::endl;

    //tri::io::Exporter<Mesh>::Save(shell, "shell_init.obj", tri::io::Mask::IOM_ALL);
}

bool ParameterizerObject::IsInitialized()
{
    return (state != ParameterizerState::InitializationFailed);
}

void ParameterizerObject::SyncChart()
{
    for (std::size_t i = 0; i < chart->fpVec.size(); ++i) {
        for (int k = 0; k < 3; ++k) {
            chart->fpVec[i]->WT(k).P() = shell.face[i].V(k)->T().P();
        }
    }
}

void ParameterizerObject::SetGradientNormTolerance(double tol)
{
    gradientNormTolerance = tol;
}

void ParameterizerObject::SetEnergyDiffTolerance(double tol)
{
    energyDiffTolerance = tol;
}

double ParameterizerObject::GetGradientNormTolerance()
{
    return gradientNormTolerance;
}

double ParameterizerObject::GetEnergyDiffTolerance()
{
    return energyDiffTolerance;
}

void ParameterizerObject::MapEnergyToShellFaceColor()
{
    energy->MapToFaceQuality(false);
    //tri::UpdateColor<Mesh>::PerFaceQualityRamp(shell, 0, 10.0 * energy->SurfaceArea());
    tri::UpdateColor<Mesh>::PerFaceQualityRamp(shell);
}

void ParameterizerObject::MapEnergyGradientToShellVertexColor()
{
    Eigen::MatrixXd G = energy->Grad();
    double maxn = 0;
    for (int i = 0; i < G.rows(); ++i) {
        maxn = std::max(maxn, G.row(i).norm());
    }
    for (int i = 0; i < G.rows(); ++i) {
        vcg::Point2d vec{G.row(i)[0], G.row(i)[1]};
        vec /= maxn;
        double angle = VecAngle(vcg::Point2d{0, 1}, vec);
        double norm = vec.Norm();
        shell.vert[i].C().SetHSVColor(angle / (2.0 * M_PI), std::pow(norm, 0.1), 1.0);
    }
}

void ParameterizerObject::MapDescentDirectionToShellVertexColor()
{
    Eigen::MatrixXd D;
    ensure_condition(opt->ComputeDescentDirection(D));
    double maxn = 0;
    for (int i = 0; i < D.rows(); ++i) {
        maxn = std::max(maxn, D.row(i).norm());
    }
    for (int i = 0; i < D.rows(); ++i) {
        vcg::Point2d vec{D.row(i)[0], D.row(i)[1]};
        vec /= maxn;
        double angle = VecAngle(vcg::Point2d{0, 1}, vec);
        double norm = vec.Norm();
        shell.vert[i].C().SetHSVColor(angle / (2.0 * M_PI), std::pow(norm, 0.25), 1.0);
    }
}

void ParameterizerObject::MapLocalGradientVarianceToShellVertexColor()
{
    /*
    Eigen::MatrixXd G = energy->Grad();
    std::vector<double> gvar(shell.VN(), 0.0);
    std::vector<double> deg(shell.VN(), 0.0);
    Eigen::Vector2d gf[3];
    for (auto& f : shell.face) {
        energy->Grad(tri::Index(shell, f), gf[0], gf[1], gf[2]);
        for (int i = 0; i < 3; ++i) {
            int vi = tri::Index(shell, f.V(i));
            Eigen::Vector2d diff = Eigen::Vector2d(G.row(vi)) - gf[i];
            gvar[vi] += diff.norm();
            deg[vi] += 1.0;
        }
    }

    double maxVar = 0.0;
    for (int i = 0; i < shell.VN(); ++i) {
        gvar[i] = gvar[i] / deg[i];
        if (gvar[i] > maxVar) maxVar = gvar[i];
    }

    for (int i = 0; i < shell.VN(); ++i) {
        double v = gvar[i] / maxVar;
        shell.vert[i].C() = vcg::Color4b(255.0, (1.0-v)*255.0, (1.0-v)*255.0, 255.0);
    }
    */
}

void ParameterizerObject::MapConformalScalingFactorsToShellVertexColor()
{
    Eigen::VectorXd csf;
    std::vector<int> coneIndices = {};
    ComputeConformalScalingFactors(csf, coneIndices);
    for (int i = 0; i < shell.VN(); ++i) {
        shell.vert[i].Q() = csf[i];
    }
    tri::UpdateColor<Mesh>::PerVertexQualityRamp(shell);
}

void ParameterizerObject::ClearShellFaceColor()
{
    for (auto& sf : shell.face) {
        sf.C() = vcg::Color4b::White;
    }
}

Mesh& ParameterizerObject::Shell()
{
    return shell;
}

ChartHandle ParameterizerObject::GetChart()
{
    return chart;
}

int ParameterizerObject::IterationCount()
{
    return iterationCount;
}

void ParameterizerObject::InitializeOptimizer()
{
    auto energy_sd = std::make_shared<SymmetricDirichletEnergy>(shell);
    switch(strategy.descent) {
    case DescentType::Gradient:
        opt = std::make_shared<GradientDescent>(energy_sd);
        break;
    case DescentType::LimitedMemoryBFGS:
        opt = std::make_shared<LBFGS>(energy_sd, 10);
        break;
    case DescentType::ScalableLocallyInjectiveMappings:
        opt = std::make_shared<SLIM>(energy_sd);
        break;
    case DescentType::CompositeMajorization:
        opt = std::make_shared<CompMaj>(energy_sd);
        break;
    default:
        ensure_condition(0 && "Unknown descent algorithm");
    }
    energy = energy_sd;
}

bool ParameterizerObject::OptimizerIsInitialized()
{
    return energy && opt;
}


/* Parameterization-related methods
 * ================================ */

bool ParameterizerObject::InitializeSolution()
{
    bool solved = false;

    if (!solved) {
        UniformSolver<Mesh> solver(shell);
        solved = solver.Solve();
    }
    if (!solved) {
        UniformSolver<Mesh> solver(shell);
        solver.SetBoundaryMapProportional();
        solved = solver.Solve();
    }

    if (!solved) {
        MeanValueSolver<Mesh> mvs(shell);
        solved = mvs.Solve();
    }

    if (!solved) {
        MeanValueSolver<Mesh> mvs(shell);
        mvs.UseCotangentWeights();
        solved = mvs.Solve();
    }

    if (solved) {
        // scale the parameterization
        double uvArea = M_PI * 0.5 * 0.5;
        double scale = std::sqrt(targetArea / uvArea);
        for (auto& sv : shell.vert) {
            sv.T().P() *= scale;
        }
        SyncShellWithUV(shell);
        state = ParameterizerState::OK;
    } else {
        std::cout << "WARNING: failed to initialize injective solution" << std::endl;
        tri::io::Exporter<Mesh>::Save(shell, "error.obj", tri::io::Mask::IOM_ALL);
        state = ParameterizerState::InitializationFailed;
    }

    return solved;
}

void ParameterizerObject::ForceWarmStart()
{
    strategy.warmStart = true;
    Reset();
}

#if 0
bool ParameterizerObject::Parameterize()
{
    if (!OptimizerIsInitialized())
        InitializeOptimizer();

    tri::io::Exporter<Mesh>::Save(shell, "shell_init.obj", tri::io::Mask::IOM_ALL);

    // check if remeshing is required during iterations
    needsRemeshing = false;
    for (auto& sf : shell.face) {
        if (sf.holeFilling) {
            needsRemeshing = true;
            break;
        }
    }

    // parameterizer state
    Timer t;
    int i;
    IterationInfo info;
    double surfaceArea = energy->SurfaceArea();
    double realSurfaceArea = energy->SurfaceAreaNotHoleFilling();
    std::vector<double> vRealEnergyVal;
    std::vector<double> vRealEnergyDiff;

    double cutTrigger = 1.25;
    double energyDiffCutTrigger = 0.02;

    vRealEnergyVal.push_back(INFINITY); // placeholder for the first iteration
    double energyCutThreshold = energy->NormalizedMinValue() * cutTrigger * realSurfaceArea;
    double energyDiffCutTriggerThreshold = energy->NormalizedMinValue() * energyDiffCutTrigger * realSurfaceArea;
    for (i = 0; i < strategy.optimizerIterations; ++i) {
        info = Iterate();

        if (info.gradientNorm < gradientNormTolerance) {
            std::cout << "Stopping because gradient magnitude is small enough (" << info.gradientNorm << ")" << std::endl;
            break;
        }

        if (info.energyDiff < (energyDiffTolerance * surfaceArea)) {
            std::cout << "Stopping because energy improvement is too small (" << info.energyDiff << ")" << std::endl;
            break;
        }

        double realEnergyVal = energy->E_IgnoreMarkedFaces(false);
        double realEnergyDiff = vRealEnergyVal.back() - realEnergyVal;
        vRealEnergyVal.push_back(realEnergyVal);
        vRealEnergyDiff.push_back(realEnergyDiff);

        /* Evaluate whether cutting is required. A cut is placed if for the last
         * 3 iterations the energy improvement was below a given threshold (which
         * hints at convergence), and at the same time the energy is far from the
         * lower bound of MIN_ENERGY * SURFACE_AREA. If this is the case, place
         * a cut and keep iterating */

        int iterationsAvailable = strategy.optimizerIterations - i;
        bool placedCut = false;
        if (strategy.applyCut && i > 20 && iterationsAvailable > 20) {
            bool guessNearConvergence = true;
            auto bit = --vRealEnergyDiff.end();
            for (int k = 0; k < 3; ++k) {
                guessNearConvergence &= *(bit--) < energyDiffCutTriggerThreshold;
            }
            if (guessNearConvergence && realEnergyVal > energyCutThreshold) {
                // Place cut
                bool success = PlaceCutWithConeSingularities(1);
                if (success) {
                    std::cout << "Placed cut with 1 cone singularity" << std::endl;
                    InitializeSolution();
                    placedCut = true;
                    vRealEnergyVal.push_back(INFINITY); // this essentially resets the cut-placing heuristic state
                } else {
                    std::cout << "Attempt to place a cut failed" << std::endl;
                }
            }
        }

        if (!placedCut && needsRemeshing && (i > 0) && (i % 30) == 0)
            RemeshHolefillingAreas();

    }
    std::cout << "Stopped after " << i << " iterations, gradient magnitude = " << info.gradientNorm
              << ", energy value = " << energy->E_IgnoreMarkedFaces() << std::endl;

    double minq, maxq;
    tri::Stat<Mesh>::ComputePerFaceQualityMinMax(shell, minq, maxq);
    std::cout << "Min distortion value = " << minq << std::endl;
    std::cout << "Max distortion value = " << maxq << std::endl;

    std::cout << "Optimization took " << t.TimeSinceLastCheck() << " seconds" << std::endl;

    return true;
}
#endif


bool ParameterizerObject::Parameterize()
{
    constexpr double CONFORMAL_SCALING_THRESHOLD = 1.98;

    if (state != ParameterizerState::OK) // bail out
        return false;

    if (!OptimizerIsInitialized())
        InitializeOptimizer();

    if (strategy.applyCut) {
        int nc = PlaceCutWithConesUntilThreshold(CONFORMAL_SCALING_THRESHOLD);
        std::cout << "Placed " << nc << " cuts" << std::endl;
    }

    if (state != ParameterizerState::OK)
        return false;

    // check if remeshing is required during iterations
    needsRemeshing = false;
    for (auto& sf : shell.face) {
        if (sf.IsHoleFilling()) {
            needsRemeshing = true;
            break;
        }
    }

    // parameterizer state
    Timer t;
    int i;
    IterationInfo info;
    double meshSurfaceArea = energy->SurfaceAreaNotHoleFilling();
    double meshCurrentEnergy = energy->E_IgnoreMarkedFaces();
    for (i = 0; i < strategy.optimizerIterations; ++i) {
        ensure_condition(state == ParameterizerState::OK);

        info = Iterate();
        if (state == ParameterizerState::IterationFailed) {
            std::cout << "Warning: failed to compute descent direction" << std::endl;
            return false;
        }

        if (info.gradientNorm < gradientNormTolerance) {
            std::cout << "Stopping because gradient magnitude is small enough (" << info.gradientNorm << ")" << std::endl;
            break;
        }
        double meshIterationEnergy = energy->E_IgnoreMarkedFaces();
        if (meshIterationEnergy >= meshCurrentEnergy) {
            std::cout << "WARNING: Stopping because energy increased" << std::endl;
            break;
        }
        double meshEnergyDiff = (meshCurrentEnergy - meshIterationEnergy);
        if (meshEnergyDiff < (energyDiffTolerance * meshSurfaceArea)) {
            std::cout << "Stopping because energy improvement is too small (" << info.energyDiff << ")" << std::endl;
            break;
        }
        meshCurrentEnergy = meshIterationEnergy;
        if (needsRemeshing && (i > 0) && (i % 30) == 0) {
            RemeshHolefillingAreas();
            opt->UpdateCache();
        }
    }
    std::cout << "Stopped after " << i << " iterations, gradient magnitude = " << info.gradientNorm
              << ", energy value = " << energy->E_IgnoreMarkedFaces() << std::endl;

    std::cout << "Optimization took " << t.TimeSinceLastCheck() << " seconds" << std::endl;

    return true;
}

IterationInfo ParameterizerObject::Iterate()
{
    if (!OptimizerIsInitialized())
        InitializeOptimizer();

    IterationInfo info = {};

    bool ok = opt->Iterate(info.gradientNorm, info.energyDiff, info.energyVal);

    if (!ok) {
        std::cout << "WARNING: Iteration failed" << std::endl;
        state = ParameterizerState::IterationFailed;
    } else {
        SyncShellWithUV(shell);
        iterationCount++;

        bool meshChanged = false;
        if (strategy.padBoundaries) {
            RemeshHolefillingAreas();
            meshChanged = true;
        }

        if (strategy.scaffold) {
            RebuildScaffold(shell, strategy.geometry, baseMesh);
            meshChanged = true;
        }

        if (meshChanged)
            opt->UpdateCache();
    }

    return info;
}

void ParameterizerObject::RemeshHolefillingAreas()
{
    ClearHoleFillingFaces(shell, true, false);
    CloseShellHoles(shell, strategy.geometry, baseMesh);
    MarkInitialSeamsAsFaux(shell, baseMesh);
}

void ParameterizerObject::PlaceCut()
{
    /*
    ComputeDistanceFromBorderOnSeams(shell);

    double maxDistance = 0;
    for (auto& sv : shell.vert) {
        if (sv.Q() < 0) ensure_condition(0);
        if (sv.Q() < INFINITY && sv.Q() > maxDistance) maxDistance = sv.Q();
    }

    MapEnergyToShellFaceColor();

    // Select candidate crease edge for cutting
    double maxEnergy = 0;
    Mesh::FacePointer startFace = nullptr;
    int cutEdge = -1;
    for (auto& sf : shell.face) {
        for (int i = 0; i < 3; ++i) {
            if (sf.IsF(i)) {
                double w = std::max(sf.V0(i)->Q(), sf.V1(i)->Q()) / maxDistance;

                if (w < INFINITY) w = 1.0;

                double weightedEnergy = energy->E(sf) * w;
                // w is INFINITY if the seam does not reach the mesh boundary
                if (std::isfinite(weightedEnergy) && weightedEnergy > maxEnergy) {
                    maxEnergy = weightedEnergy;
                    startFace = &sf;
                    cutEdge = i;
                }
            }
        }
    }
    ensure_condition(startFace != nullptr);

    //tri::io::Exporter<Mesh>::Save(shell, "shell_before_cut.obj", tri::io::Mask::IOM_ALL);

    std::cout << "Selected face " << tri::Index(shell, startFace) << " edge " << cutEdge << std::endl;

    PosF p(startFace, cutEdge);
    ensure_condition(!p.IsBorder());

    tri::UpdateFlags<Mesh>::FaceClearFaceEdgeS(shell);
    PosF boundaryPos = SelectShortestSeamPathToBoundary(shell, p);
    //SelectShortestSeamPathToPeak(shell, p);
    bool corrected = RectifyCut(shell, boundaryPos);

    tri::CutMeshAlongSelectedFaceEdges(shell);
    CleanupShell(shell);

    tri::UpdateTopology<Mesh>::FaceFace(shell);

    opt->UpdateCache();
    //tri::io::Exporter<Mesh>::Save(shell, "shell_after_cut.obj", tri::io::Mask::IOM_ALL);
    */
}

/* What this method does is:
 *  1. Remove any hole filling faces from the shell
 *  2. Sync the shell with the model
 *  3. Close the holes
 *  4. Compute 'ncones' cone singularities using the method described in
 *     Ben-Chen et al. 'Conformal flattening by curvature prescription and
 *     metric scaling' (2008), restricted to the vertices that lie on uv seams
 *  5. Iteratively cut the shell, connecting each cone with the boundary
 *  6. Remove the hole filling faces, sync the shell with its UV configuration
 *     and close the holes in UV space */
bool ParameterizerObject::PlaceCutWithConeSingularities(int ncones)
{
    ensure_condition(ncones > 0);

    // sanity check
    ComputeDistanceFromBorderOnSeams(shell);
    bool seamsToBoundary = false;
    for (auto& sv : shell.vert) {
        if (!sv.IsB() && sv.Q() < INFINITY) {
            seamsToBoundary = true;
            break;
        }
    }
    if (!seamsToBoundary)
        return false;

    // Put the shell in a configuration suitable to compute the conformal scaling factors
    ClearHoleFillingFaces(shell, strategy.padBoundaries, strategy.scaffold);

    SyncShellWithModel(shell, baseMesh);

    ComputeBoundaryInfo(shell);
    CloseMeshHoles(shell);
    MarkInitialSeamsAsFaux(shell, baseMesh);

    std::vector<int> coneIndices;
    FindCones(ncones, coneIndices);

    ComputeDistanceFromBorderOnSeams(shell);

    // At this point the cones are guaranteed to be placed on a mesh vertex
    tri::UpdateFlags<Mesh>::VertexClearS(shell);
    for (int i : coneIndices)
        shell.vert[i].SetS();

    std::vector<PosF> pv;
    for (auto& sf : shell.face) {
        if (sf.IsMesh()) {
            for (int i = 0; i < 3; ++i) {
                if (sf.IsF(i) && (sf.V0(i)->IsS() || sf.V1(i)->IsS())) {
                    PosF p(&sf, i);
                    p.V()->ClearS();
                    p.VFlip()->ClearS();
                    pv.push_back(p);
                }
            }
        }
    }

    for (auto p : pv) {
        std::cout << "Selected face " << tri::Index(shell, p.F()) << " edge " << p.E() << std::endl;
        tri::UpdateFlags<Mesh>::FaceClearFaceEdgeS(shell);
        PosF boundaryPos = SelectShortestSeamPathToBoundary(shell, p);
        //SelectShortestSeamPathToPeak(shell, p);
        bool corrected = RectifyCut(shell, boundaryPos);
        tri::CutMeshAlongSelectedFaceEdges(shell);
        CleanupShell(shell);
        tri::UpdateTopology<Mesh>::FaceFace(shell);

        ensure_condition(tri::Clean<Mesh>::CountConnectedComponents(shell) == 1);
    }

    ClearHoleFillingFaces(shell, true, false);
    ComputeBoundaryInfo(shell); // boundary changed after the cut
    SyncShellWithUV(shell);

    if (strategy.padBoundaries)
        CloseShellHoles(shell, strategy.geometry, baseMesh);

    if (OptimizerIsInitialized())
        opt->UpdateCache();

    return true;
}

#include <wrap/io_trimesh/export.h>

int ParameterizerObject::PlaceCutWithConesUntilThreshold(double conformalScalingThreshold)
{
    if (state != ParameterizerState::OK)
        return 0;
    /*
    // sanity check
    if (strategy.scaffold) {
        ClearHoleFillingFaces(shell, false, true);
    }
    ComputeDistanceFromBorderOnSeams(shell);
    bool seamsToBoundary = false;
    for (auto& sv : shell.vert) {
        if (!sv.IsB() && sv.Q() < INFINITY) {
            seamsToBoundary = true;
            break;
        }
    }
    if (!seamsToBoundary) {
        if (strategy.scaffold) {
            BuildScaffold(shell, strategy.geometry, baseMesh);
            if (OptimizerIsInitialized())
                opt->UpdateCache();

        }
        return 0;
    }
    */

    static int cut_shell = 0;

    // Put the shell in a configuration suitable to compute the conformal scaling factors
    ClearHoleFillingFaces(shell, true, true);

    std::vector<double> savedUVs;
    for (auto& sf : shell.face) {
        for (int i = 0; i < 3; ++i) {
            savedUVs.push_back(sf.P(i).X());
            savedUVs.push_back(sf.P(i).Y());
        }
    }

    SyncShellWithModel(shell, baseMesh);
    for (auto& sf : shell.face) {
        for (int i = 0; i < 3; ++i) {
            if (sf.IsF(i)) {
                sf.C() = vcg::Color4b::Blue;
                sf.FFp(i)->C() = vcg::Color4b::Blue;
            }
        }
    }

    tri::Clean<Mesh>::RemoveDuplicateVertex(shell);
    tri::UpdateTopology<Mesh>::FaceFace(shell);
    int splitCount;
    while ((splitCount = tri::Clean<Mesh>::SplitNonManifoldVertex(shell, 0.15)) > 0) {
        tri::Allocator<Mesh>::CompactEveryVector(shell);
    }
    tri::Allocator<Mesh>::CompactEveryVector(shell);

    ComputeBoundaryInfo(shell);
    CloseMeshHoles(shell);
    MarkInitialSeamsAsFaux(shell, baseMesh);

    int numPlacedCones = 0;

    static    int cut_id = 0;
    for (;;) {
        std::vector<int> coneIndices;
        FindConesWithThreshold(conformalScalingThreshold, coneIndices);

        if (coneIndices.size() > 0) {
            // Cone singularities were found, place the cut
            // The cones are guaranteed to be placed on a mesh vertex since
            // CloseMeshHoles() fills holes by closing ears (no vertices are added)

            ComputeDistanceFromBorderOnSeams(shell);
            // Select cone vertices
            tri::UpdateFlags<Mesh>::VertexClearS(shell);
            for (int i : coneIndices)
                shell.vert[i].SetS();
            // Store Pos objects that lie on seams and point to cone vertices (1 pos per cone)
            std::vector<PosF> pv;
            for (auto& sf : shell.face) {
                if (sf.IsMesh()) {
                    for (int i = 0; i < 3; ++i) {
                        if (sf.IsF(i) && (sf.V0(i)->IsS() || sf.V1(i)->IsS())) {
                            PosF p(&sf, i);
                            p.V()->ClearS();
                            p.VFlip()->ClearS();
                            pv.push_back(p);
                        }
                    }
                }
            }
            // Repeatedly cut starting from each pos
            for (auto p : pv) {
                std::cout << "Selected face " << tri::Index(shell, p.F()) << " edge " << p.E() << std::endl;
                tri::UpdateFlags<Mesh>::FaceClearFaceEdgeS(shell);
                PosF boundaryPos = SelectShortestSeamPathToBoundary(shell, p);
                //SelectShortestSeamPathToPeak(shell, p);
                bool corrected = RectifyCut(shell, boundaryPos);

                /*
                {
                    MyMesh em;

                    for (auto& sf : shell.face) {
                        for (int i = 0; i < 3; ++i) {
                            if (sf.IsFaceEdgeS(i)) {
                                tri::Allocator<MyMesh>::AddEdge(em, sf.P0(i), sf.P1(i));
                            }
                        }
                    }
                    tri::Clean<MyMesh>::RemoveDuplicateEdge(em);
                    tri::Clean<MyMesh>::RemoveDuplicateVertex(em);
                    tri::Allocator<MyMesh>::CompactEveryVector(em);

                    stringstream ss;
                    ss << "cut_" << cut_id++ << ".obj";
                    tri::io::Exporter<MyMesh>::Save(em, ss.str().c_str());
                }
                */

                tri::CutMeshAlongSelectedFaceEdges(shell);
                CleanupShell(shell);
                tri::UpdateTopology<Mesh>::FaceFace(shell);
            }

            numPlacedCones += (int) coneIndices.size();
        } else {
            break;
        }
    }


    // Cleanup (and restore hole-filling and scaffold faces if needed)
    ClearHoleFillingFaces(shell, true, true);

    /*
    auto ia = GetFaceIndexAttribute(shell);
    for (auto& sf : shell.face) {
        assert(ia[sf] != -1);
        auto& mf = baseMesh.face[ia[sf]];
        for (int i = 0; i < 3; ++i) {
            sf.WT(i) = mf.WT(i);
        }
    }
    tri::io::Exporter<Mesh>::Save(shell, "shell.obj", tri::io::Mask::IOM_WEDGTEXCOORD);
    */

    // restore the vertex texture coordinates
    ensure_condition(shell.FN() == ((int) savedUVs.size() / 6));
    double *coordPtr = savedUVs.data();
    for (auto& sf : shell.face) {
        for (int i = 0; i < 3; ++i) {
            sf.V(i)->T().P().X() = *coordPtr++;
            sf.V(i)->T().P().Y() = *coordPtr++;
        }
    }

    ComputeBoundaryInfo(shell); // boundary changed after the cut
    SyncShellWithUV(shell);

    if (strategy.padBoundaries) {
        CloseShellHoles(shell, strategy.geometry, baseMesh);
    }


    // Initialize solution if cones were placed
    if (numPlacedCones > 0)
        InitializeSolution();

    if (strategy.scaffold)
        BuildScaffold(shell, strategy.geometry, baseMesh);

    MarkInitialSeamsAsFaux(shell, baseMesh);

    if (OptimizerIsInitialized())
        opt->UpdateCache();

    return numPlacedCones;
}

void ParameterizerObject::FindCones(int ncones, std::vector<int>& coneIndices)
{
    coneIndices.clear();

    ComputeDistanceFromBorderOnSeams(shell);

    for (int i = 0; i < ncones; ++i) {
        Eigen::VectorXd csf;
        if (!ComputeConformalScalingFactors(csf, coneIndices)) {
            std::cout << "Warning: scaling factors computation failed" << std::endl;
            return;
        }

        // Select seam vertices
        tri::UpdateFlags<Mesh>::VertexClearS(shell);
        for (auto& sv : shell.vert)
            if (sv.Q() < INFINITY)
                sv.SetS();

        // Find the max scale factor at seam vertices
        double maxscale = 0.0;
        double minscale = std::numeric_limits<double>::infinity();
        int max_k = -1;
        int min_k = -1;
        for (int k = 0; k < shell.VN(); ++k) {
            if (shell.vert[k].IsS()) {
                double s = std::abs(csf[k]);
                if (s > maxscale) {
                    maxscale = s;
                    max_k = k;
                }
                if (s < minscale) {
                    minscale = s;
                    min_k = k;
                }
            }
        }

        std::cout << "CONFORMAL SCALING FACTORS: max = " << maxscale << " | min = " << minscale <<  " | ratio = " << minscale / maxscale << std::endl;

        coneIndices.push_back(max_k);
    }
}

/* This function is a bit messed up, in theory one could choose an arbitrary number
 * of cone vertices, but it makes little sense to me to treat a 'designated' cone
 * as boundary vertex without actually cutting the surface. So I restrict the number
 * of cones to 1, forcing to cut the surface repeatedly. */
void ParameterizerObject::FindConesWithThreshold(double conformalScalingThreshold, std::vector<int>& coneIndices)
{
    constexpr int MAX_NUM_CONES = 1;

    coneIndices.clear();

    ComputeDistanceFromBorderOnSeams(shell);
    bool seamsToBoundary = false;
    for (auto& sv : shell.vert) {
        if (!sv.IsB() && sv.Q() < INFINITY) {
            seamsToBoundary = true;
            break;
        }
    }

    // bail out if no seams reach the mesh boundary
    if (!seamsToBoundary)
        return;

    for (int i = 0; i < MAX_NUM_CONES; ++i) {
        Eigen::VectorXd csf;
        if(!ComputeConformalScalingFactors(csf, coneIndices)) {
            std::cout << "Warning: scaling factors computation failed" << std::endl;
            return;
        }

        // Select seam vertices that reach the boundary
        tri::UpdateFlags<Mesh>::VertexClearS(shell);
        for (auto& sv : shell.vert)
            if (sv.Q() > 0 && sv.Q() < INFINITY)
                sv.SetS();

        // Find the max scale factor at seam vertices
        double maxscale = 0.0;
        double minscale = std::numeric_limits<double>::infinity();
        int max_k = -1;
        int min_k = -1;
        for (int k = 0; k < shell.VN(); ++k) {
            if (shell.vert[k].IsS()) {
                double s = std::abs(csf[k]);
                if (s > maxscale) {
                    maxscale = s;
                    max_k = k;
                }
                if (s < minscale) {
                    minscale = s;
                    min_k = k;
                }
            }
        }

        std::cout << "CONFORMAL SCALING FACTORS: max = " << maxscale << " | min = " << minscale <<  " | ratio = " << minscale / maxscale << std::endl;

        if (maxscale > conformalScalingThreshold)
            coneIndices.push_back(max_k);
        else
            break;
    }

}

/* assumes that cone singularity vertices have the selected flag set */
static void ComputeMarkovMatrix(Mesh& m, Eigen::SparseMatrix<double>& L, SimpleTempData<Mesh::VertContainer, int>& Idx)
{
    using Td = Eigen::Triplet<double>;

    L.resize(m.VN(), m.VN());
    L.setZero();
    std::unordered_map<Mesh::VertexPointer, std::unordered_set<Mesh::VertexPointer>> vadj;
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            vadj[f.V0(i)].insert(f.V1(i));
            vadj[f.V0(i)].insert(f.V2(i));
        }
    }
    std::vector<Td> tri;
    for (auto& entry: vadj) {
        if (!entry.first->IsS()) {
            double c = entry.second.size();
            for (auto & vp : entry.second) {
                tri.push_back(Td(Idx[entry.first], Idx[vp], (1.0 / c)));
            }
        } else {
            tri.push_back(Td(Idx[entry.first], Idx[entry.first], 1.0));
        }
    }
    L.setFromTriplets(tri.begin(), tri.end());
    L.makeCompressed();
}

void ComputeCotangentWeightedLaplacian(Mesh& shell, Eigen::SparseMatrix<double>& L, SimpleTempData<Mesh::VertContainer,int>& Idx)
{
    using Td = Eigen::Triplet<double>;

    L.resize(shell.VN(), shell.VN());
    L.setZero();
    std::vector<Td> tri;
    double eps = std::numeric_limits<double>::epsilon();
    for (auto &sf : shell.face) {
        for (int i = 0; i < 3; ++i) {
            Mesh::VertexPointer vi = sf.V0(i);
            Mesh::VertexPointer vj = sf.V1(i);
            Mesh::VertexPointer vk = sf.V2(i);

            ensure_condition(Idx[vi] >= 0 && Idx[vi] < shell.VN());
            ensure_condition(Idx[vj] >= 0 && Idx[vj] < shell.VN());
            ensure_condition(Idx[vk] >= 0 && Idx[vk] < shell.VN());

            double alpha_ij = std::max(VecAngle(vi->P() - vk->P(), vj->P() - vk->P()), eps); // angle at k
            double alpha_ik = std::max(VecAngle(vi->P() - vj->P(), vk->P() - vj->P()), eps); // angle at j
            double weight_ij = 0.5 * std::tan(M_PI_2 - alpha_ij);
            double weight_ik = 0.5 * std::tan(M_PI_2 - alpha_ik);

            if (!std::isfinite(weight_ij))
                weight_ij = 1e-8;

            if (!std::isfinite(weight_ik))
                weight_ik = 1e-8;

            ensure_condition(std::isfinite(weight_ij));
            ensure_condition(std::isfinite(weight_ik));

            tri.push_back(Td(Idx[vi], Idx[vj], -weight_ij));
            tri.push_back(Td(Idx[vi], Idx[vk], -weight_ik));
            tri.push_back(Td(Idx[vi], Idx[vi], (weight_ij + weight_ik)));
        }
    }
    L.setFromTriplets(tri.begin(), tri.end());
    L.makeCompressed();
}


/* Reference: BenChen 2009 curvature prescription metric scaling */
bool ParameterizerObject::ComputeConformalScalingFactors(Eigen::VectorXd& csf, const std::vector<int>& coneIndices)
{
    Mesh& m = shell;

    Timer timer;

    std::cout << "Computing conformal scaling factors..." << std::endl;

    // Permutate indices so that boundary vertices are the last ones
    SimpleTempData<Mesh::VertContainer, int> index{m.vert};
    tri::UpdateFlags<Mesh>::VertexBorderFromFaceAdj(m);
    tri::UpdateFlags<Mesh>::VertexClearS(m);
    int n_internal = 0;
    int j = m.VN() - 1;

    for (auto i : coneIndices) {
        auto& v = m.vert[i];
        v.SetS();
        index[v] = j--;
    }

    for (auto& v : m.vert) {
        if (!v.IsS()) {
            if (v.IsB()) {
                index[v] = j--;
                v.SetS();
            } else
                index[v] = n_internal++;
        }
    }
    int n_boundary = m.VN() - n_internal;

    // Compute Laplacian and Markov system
    Eigen::SparseMatrix<double> M;

    Eigen::VectorXd Korig{m.VN()};
    Korig.setZero();
    for (auto& sf : m.face) {
        //auto& f = baseMesh.face[ia[sf]];
        for (int i = 0; i < 3; ++i) {
            double angle = VecAngle(sf.P1(i) - sf.P0(i), sf.P2(i) - sf.P0(i));
            if (std::isfinite(angle))
                Korig[index[sf.V(i)]] -= angle;
        }
    }
    Korig.segment(0, n_internal) += Eigen::VectorXd::Constant(n_internal, 2.0 * M_PI);
    Korig.segment(n_internal, n_boundary) += Eigen::VectorXd::Constant(n_boundary, M_PI);

    //std::cout << "Sum Korig = " << Korig.sum() << std::endl;

    // Compute the metric scaling of the new target curvature
    ComputeCotangentWeightedLaplacian(m, M, index);
    M += Eigen::VectorXd::Constant(M.rows(), 1e-12).asDiagonal();

    Eigen::SparseMatrix<double> M_internal = M.block(0, 0, n_internal, n_internal);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverLDLT;
    solverLDLT.compute(M_internal);
    if (solverLDLT.info() != Eigen::Success) {
        std::cout << "WARNING ComputeConformalScalingFactors(): preliminary factorization failed" << std::endl;
        return false;
    }

    // Compute the new curvature vector Knew
    Eigen::VectorXd Knew = Eigen::VectorXd::Zero(m.VN());
    Eigen::VectorXd G(n_internal);

    // This loops takes a lot of time (a backward solve for each boundary vertex...)
    for (int i = 0; i < n_boundary; ++i) {
        G = solverLDLT.solve(-M.block(0, n_internal+i, n_internal, 1));
        if (solverLDLT.info() != Eigen::Success) {
            std::cout << "WARNING ComputeConformalScalingFactors(): backward substitution to compute G(" << i << ") failed" << std::endl;
            return false;
        }
        Knew[n_internal+i] = Korig[n_internal+i] + G.dot(Korig.head(n_internal));
    }

    //std::cout << "Sum Knew = " << Knew.sum() << std::endl;
    std::cout << "Computation of the target curvatures took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;

    solverLDLT.compute(M);
    if (solverLDLT.info() != Eigen::Success) {
        std::cout << "WARNING ComputeConformalScalingFactors(): factorization of the cotangent-weighted Laplacian failed" << std::endl;
        return false;
    }

    Eigen::VectorXd p = Knew - Korig;

    Eigen::VectorXd scalingFactors = solverLDLT.solve(Knew - Korig);
    csf.resize(scalingFactors.size());
    for (int i = 0; i < m.VN(); ++i) {
        csf[i] = scalingFactors[index[m.vert[i]]];
    }

    double maxscale = 0.0;
    double minscale = std::numeric_limits<double>::infinity();
    for (int i = 0; i < m.VN(); ++i) {
        double s = std::abs(csf[i]);
        //ensure_condition(scalingFactors[i] > 0.0);
        if (s > maxscale) {
            maxscale = s;
        }
        if (s < minscale) {
            minscale = s;
        }
    }

    std::cout << "Solution of the scaling factors system took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    std::cout << "Overall time: " << timer.TimeElapsed() << " seconds" << std::endl;

    return true;
}


