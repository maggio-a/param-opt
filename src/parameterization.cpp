#include "parameterization.h"
#include "mesh_graph.h"
#include "energy.h"
#include "iterative_solvers.h"
#include "mesh_utils.h"
#include "math_utils.h"
#include "timer.h"
#include "linear_solvers.h"
#include "uv.h"
#include "polygon2_triangulator.h"
#include "logging.h"

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

/* ! FIXME useExistingUV is ignored, remove ! */
bool BuildShell(Mesh& shell, FaceGroup& fg, ParameterizationGeometry targetGeometry, bool useExistingUV)
{
    static double t0 = 0;
    static double t1 = 0;
    static double t2 = 0;
    static double t3 = 0;
    static double t4 = 0;
    static double t5 = 0;
    static double t6 = 0;

    Timer t;

    Mesh& m = fg.mesh;

    CopyToMesh(fg, shell);

    t0 += t.TimeElapsed();

    tri::Clean<Mesh>::RemoveDuplicateVertex(shell);
    tri::UpdateBounding<Mesh>::Box(shell);
    tri::UpdateTopology<Mesh>::FaceFace(shell);

    t1 += t.TimeElapsed();

    int splitCount;
    while ((splitCount = tri::Clean<Mesh>::SplitNonManifoldVertex(shell, 0.15)) > 0)
        ;
    tri::Allocator<Mesh>::CompactEveryVector(shell);

    t2 += t.TimeElapsed();

    // First attempt, topological cut
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
    }

    t3 += t.TimeElapsed();

    // Second attempt, pick a face and cut along two edges
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

    t4 += t.TimeElapsed();

    // Failed to produce a parameterizable shell
    if (!Parameterizable(shell))
        return false;

    CloseMeshHoles(shell);

    t5 += t.TimeElapsed();

    // Compute the target shapes for the shell faces
    auto psi = GetParameterizationScaleInfoAttribute(m);
    double avg3D = psi().surfaceArea / psi().numNonZero;
    double avgUV = psi().parameterArea / psi().numNonZero;

    double targetArea = 0;

    auto sa = GetShell3DShapeAttribute(shell);
    auto ia = GetFaceIndexAttribute(shell);
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
        } else if (sf.IsMesh()) {
            auto& mf = m.face[ia[sf]];
            if (targetGeometry == Model) {
                target.P[0] = mf.P(0);
                target.P[1] = mf.P(1);
                target.P[2] = mf.P(2);
            } else if (targetGeometry == Texture) {
                // Interpolate between texture and mesh face shapes to mitigate distortion
                const Point2d& u0 = wtcsa[mf].tc[0].P();
                const Point2d& u1 = wtcsa[mf].tc[1].P();
                const Point2d& u2 = wtcsa[mf].tc[2].P();
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

                    // scale using the average scaling factor rather than the one of
                    // the triangle, since heavily distorted triangle are likely to
                    // have areas that are too small
                    x10 *= psi().scale;
                    x20 *= psi().scale;

                    // compute the singular values of the transformation matrix, s2 > s1
                    // ref for the formula: Smith&Schaefer 2015
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

                    interpolationFactor = 0;

                    Point2d v10 = (1.0 - interpolationFactor) * u10 + (interpolationFactor) * x10;
                    Point2d v20 = (1.0 - interpolationFactor) * u20 + (interpolationFactor) * x20;

                    target.P[0] = Point3d(0, 0, 0);
                    target.P[1] = Point3d(v10[0], v10[1], 0);
                    target.P[2] = Point3d(v20[0], v20[1], 0);
                }
            }
        } else {
            ensure_condition(0 && "Unexpected face type when building shell");
        }
        tsa[sf] = target;
        targetArea += ((target.P[1] - target.P[0]) ^ (target.P[2] - target.P[0])).Norm() / 2.0;

        bool ok = std::isfinite(targetArea) && (targetArea > 0);
        if (!ok) {
            LOG_DEBUG << "Problem with face " << tri::Index(shell, sf);
            tri::io::Exporter<Mesh>::Save(shell, "non_finite_target_area.obj", tri::io::Mask::IOM_ALL);
        }
        ensure_condition(std::isfinite(targetArea));
        ensure_condition(targetArea > 0);

        sa[sf].P[0] = sf.P(0);
        sa[sf].P[1] = sf.P(1);
        sa[sf].P[2] = sf.P(2);
    }

    t6 += t.TimeElapsed();

    LOG_DEBUG << "\t\t\t[TIMING BuildShell()] t0 = " << t0;
    LOG_DEBUG << "\t\t\t[TIMING BuildShell()] t1 = " << t1;
    LOG_DEBUG << "\t\t\t[TIMING BuildShell()] t2 = " << t2;
    LOG_DEBUG << "\t\t\t[TIMING BuildShell()] t3 = " << t3;
    LOG_DEBUG << "\t\t\t[TIMING BuildShell()] t4 = " << t4;
    LOG_DEBUG << "\t\t\t[TIMING BuildShell()] t5 = " << t5;
    LOG_DEBUG << "\t\t\t[TIMING BuildShell()] t6 = " << t6;

    return true;
}

void SyncShellWithUV(Mesh& shell)
{
    for (auto& v : shell.vert) {
        v.P().X() = v.T().U();
        v.P().Y() = v.T().V();
        v.P().Z() = 0.0;
    }
    tri::UpdateBounding<Mesh>::Box(shell);
}

// TODO this method should CALL ClearHoleFillingFaces(shell, true, true) otherwise it's unsafe
void SyncShellWith3D(Mesh& shell)
{
    auto sa = GetShell3DShapeAttribute(shell);
    for (auto& sf : shell.face) {
        ensure_condition(sf.IsMesh());
        for (int i = 0; i < 3; ++i)
            sf.P(i) = sa[sf].P[i];
    }
    tri::UpdateBounding<Mesh>::Box(shell);
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


/* ParameterizerObject class implementation
 * ======================================== */

std::ostream& operator<<(std::ostream& os, const ParameterizerObject::Status& status)
{
    std::string str;
    switch(status) {
    case ParameterizerObject::Status::OK:
        str = "OK";
        break;
    case ParameterizerObject::Status::Error_InitFailed:
        str = "Err_InitFailed";
        break;
    case ParameterizerObject::Status::Error_IterationFailed:
        str = "Err_IterFailed";
        break;
    case ParameterizerObject::Status::Error_UnfeasibleShell:
        str = "Err_Unfeasible";
        break;
    case ParameterizerObject::Status::Initialized:
        str = "Initialized";
        break;
    case ParameterizerObject::Status::NoInit:
        str = "NoInit";
        break;
    }
    os << str;
    return os;
}

ParameterizerObject::ParameterizerObject(ChartHandle c, ParameterizationStrategy strat)
    : shell{},
      baseMesh{c->mesh},
      chart{c},
      strategy(strat),
      energy{},
      opt{},
      iterationCount{0},
      gradientNormTolerance{1e-3},
      energyDiffTolerance{1e-6},
      conformalScalingThreshold{1.98},
      targetArea{0},
      status{NoInit}
{
    shell.Clear();
    energy = nullptr;
    opt = nullptr;
    iterationCount = 0;
    gradientNormTolerance = 1e-3;
    energyDiffTolerance = 1e-6;
    conformalScalingThreshold = 1.98;
    targetArea = 0;
    SetStatus(NoInit);
}

ParameterizerObject::~ParameterizerObject()
{
    // empty
}

void ParameterizerObject::SetStatus(Status s)
{
    if (!ErrorState())
        status = s;
}

bool ParameterizerObject::ErrorState()
{
    return (status == Error_InitFailed || status == Error_IterationFailed || status == Error_UnfeasibleShell);
}

void ParameterizerObject::Initialize()
{
    if (status == Initialized)
        return;

    if (!BuildShell(shell, *chart, strategy.geometry, strategy.warmStart)) {
        SetStatus(Error_UnfeasibleShell);
        return;
    }

    MarkInitialSeamsAsFaux(shell, baseMesh);

    // cut the mesh (if strategy allows) (3D)
    if (strategy.applyCut) {
        LOG_DEBUG << "Cutting chart...";
        int nc = PlaceCutWithConesUntilThreshold_3D(conformalScalingThreshold);
        LOG_DEBUG << "Placed " << nc << " cuts";
    }

    LOG_DEBUG << "Computing initial solution";
    // initialize solution (3D)
    if (!InitializeSolution()) {
        SetStatus(Error_InitFailed);
        return;
    }

    // clear hole-filling faces (3D)
    ClearHoleFillingFaces(shell, true, true);

    // sync shell with uv (3D -> 2D)
    SyncShellWithUV(shell);

    // close shell holes (2D)
    if (strategy.padBoundaries) {
        LOG_DEBUG << "Closing UV holes";
        CloseShellHoles(shell, strategy.geometry, baseMesh);
    }

    // build scaffold (2D)
    if (strategy.scaffold) {
        LOG_DEBUG << "Building initial scaffold";
        BuildScaffold(shell, strategy.geometry, baseMesh);
    }

    auto tsa = GetTargetShapeAttribute(shell);
    targetArea = 0;
    for (auto& sf : shell.face) {
        CoordStorage target = tsa[sf];
        targetArea += ((target.P[1] - target.P[0]) ^ (target.P[2] - target.P[0])).Norm() / 2.0;
    }

    LOG_DEBUG << "Target surface area = " << targetArea;

    LOG_DEBUG << "Initializing optimizer";
    InitializeOptimizer();

    SetStatus(Initialized);
    LOG_DEBUG << "Initialization of ParameterizerObject complete";
}

ParameterizerObject::Status ParameterizerObject::GetStatus()
{
    return status;
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
        solved = (solver.Solve() && CheckLocalInjectivity(shell));
        if (!solved)
            tri::io::Exporter<Mesh>::Save(shell, "error_tutte.obj", tri::io::Mask::IOM_ALL);
    }
    if (!solved) {
        UniformSolver<Mesh> solver(shell);
        solver.SetBoundaryMapProportional();
        solved = (solver.Solve() && CheckLocalInjectivity(shell));
        if (!solved)
            tri::io::Exporter<Mesh>::Save(shell, "error_tutte_prop.obj", tri::io::Mask::IOM_ALL);
    }
    if (!solved) {
        MeanValueSolver<Mesh> mvs(shell);
        solved = (mvs.Solve() && CheckLocalInjectivity(shell));
        if (!solved)
            tri::io::Exporter<Mesh>::Save(shell, "error_mvs.obj", tri::io::Mask::IOM_ALL);
    }
    if (!solved) {
        MeanValueSolver<Mesh> mvs(shell);
        mvs.UseCotangentWeights();
        solved = (mvs.Solve() && CheckLocalInjectivity(shell));
        if (!solved)
            tri::io::Exporter<Mesh>::Save(shell, "error_cotan.obj", tri::io::Mask::IOM_ALL);
    }

    // Scale the shell size to match the target shape area
    if (solved) {
        double shellUvArea = 0;
        double shellTargetArea = 0;
        auto tsa = GetTargetShapeAttribute(shell);
        for (auto& sf : shell.face) {
            shellUvArea += ((sf.V(1)->T().P() - sf.V(0)->T().P()) ^ (sf.V(2)->T().P() - sf.V(0)->T().P())) / 2.0;
            shellTargetArea += ((tsa[sf].P[1] - tsa[sf].P[0]) ^ (tsa[sf].P[2] - tsa[sf].P[0])).Norm() / 2.0;
        }

        double areaScaleFactor = std::sqrt(shellTargetArea / shellUvArea);
        ensure_condition(areaScaleFactor > 0);

        for (auto& v : shell.vert)
            v.T().P() *= areaScaleFactor;
    }

    return solved;
}

bool ParameterizerObject::Parameterize()
{
    if (status != Initialized) {
        LOG_DEBUG << "ParameterizerObject is not initialized (" << GetStatus() << ")";
        return false;
    }

    if (!OptimizerIsInitialized())
        InitializeOptimizer();

    // check if remeshing is required during iterations
    bool needsRemeshing = false;
    for (auto& sf : shell.face) {
        if (sf.IsHoleFilling()) {
            needsRemeshing = true;
            break;
        }
    }

    SetStatus(OK);

    if (strategy.padBoundaries) {
        int nh = tri::Clean<Mesh>::CountHoles(shell);
        ensure_condition(nh == 1);
    }

    // parameterizer state
    Timer t;
    int i;
    IterationInfo info;
    double meshSurfaceArea = energy->SurfaceAreaNotHoleFilling();
    double meshCurrentEnergy = energy->E_IgnoreMarkedFaces();
    for (i = 0; i < strategy.optimizerIterations; ++i) {
        ensure_condition(status == OK);

        info = Iterate();
        if (status == Error_IterationFailed) {
            LOG_WARN << "Failed to compute descent direction";
            return false;
        }

        if (info.gradientNorm < gradientNormTolerance) {
            LOG_VERBOSE << "Stopping because gradient magnitude is small enough (" << info.gradientNorm << ")";
            break;
        }
        double meshIterationEnergy = energy->E_IgnoreMarkedFaces();
        if (meshIterationEnergy >= meshCurrentEnergy) {
            LOG_VERBOSE << "Stopping because energy increased";
            break;
        }
        double meshEnergyDiff = (meshCurrentEnergy - meshIterationEnergy);
        if (meshEnergyDiff < (energyDiffTolerance * meshSurfaceArea)) {
            LOG_VERBOSE << "Stopping because energy improvement is too small (" << info.energyDiff << ")";
            break;
        }
        meshCurrentEnergy = meshIterationEnergy;
        if (needsRemeshing && (i > 0) && (i % 30) == 0) {
            RemeshHolefillingAreas();
            opt->UpdateCache();
        }
    }
    LOG_INFO << "Stopped after " << i << " iterations, gradient magnitude = " << info.gradientNorm
             << ", normalized energy value = " << (energy->E_IgnoreMarkedFaces() / meshSurfaceArea);

    LOG_INFO << "Optimization took " << t.TimeSinceLastCheck() << " seconds";

    return true;
}

IterationInfo ParameterizerObject::Iterate()
{
    if (!OptimizerIsInitialized())
        InitializeOptimizer();

    IterationInfo info = {};

    bool ok = opt->Iterate(info.gradientNorm, info.energyDiff, info.energyVal);

    if (!ok) {
        LOG_WARN << "Iteration failed";
        SetStatus(Error_IterationFailed);
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











#include "cones.h"
int ParameterizerObject::PlaceCutWithApproximateCone()
{
    if (ErrorState())
        return 0;

    // Put the shell in a configuration suitable to compute the conformal scaling factors
    ClearHoleFillingFaces(shell, true, true);

    SyncShellWith3D(shell);
    tri::UpdateTopology<Mesh>::FaceFace(shell);

    CloseMeshHoles(shell);


    MarkInitialSeamsAsFaux(shell, baseMesh);

    {
        unsigned cone = FindApproximateCone(shell, true);

        LOG_DEBUG << "Cone at vertex " << cone;

        ComputeDistanceFromBorderOnSeams(shell);

        // Select cone vertices
        tri::UpdateFlags<Mesh>::VertexClearS(shell);
        shell.vert[cone].SetS();
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
            LOG_DEBUG << "Selected face " << tri::Index(shell, p.F()) << " edge " << p.E();
            tri::UpdateFlags<Mesh>::FaceClearFaceEdgeS(shell);
            PosF boundaryPos = SelectShortestSeamPathToBoundary(shell, p);
            //SelectShortestSeamPathToPeak(shell, p);
            RectifyCut(shell, boundaryPos);
            tri::CutMeshAlongSelectedFaceEdges(shell);
            CleanupShell(shell);
            tri::UpdateTopology<Mesh>::FaceFace(shell);
        }
    }

    MarkInitialSeamsAsFaux(shell, baseMesh);

    // need to update the boundary information, which is implicitly computed
    // by the hole filling function (ugly...)
    ClearHoleFillingFaces(shell, true, true);
    CloseMeshHoles(shell);

    int numPlacedCones = 1;

    // if cones were found we reinitialize the solution, otherwise just leave
    // everything as is. In this second case we can simply return the shell to
    // the previous uv state, which is not changed.

    if (numPlacedCones > 0) {
        InitializeSolution();
    }

    // Cleanup (and restore hole-filling and scaffold faces if needed)
    ClearHoleFillingFaces(shell, true, true);
    SyncShellWithUV(shell);

    if (strategy.padBoundaries)
        CloseShellHoles(shell, strategy.geometry, baseMesh);

    if (strategy.scaffold) {
        BuildScaffold(shell, strategy.geometry, baseMesh);
    }

    MarkInitialSeamsAsFaux(shell, baseMesh);

    if (OptimizerIsInitialized()) {
        opt->UpdateCache();
    }

    return numPlacedCones;
}
















/* This function places the cuts with the shell is in its 3D configuration */
int ParameterizerObject::PlaceCutWithConesUntilThreshold_3D(double cst)
{
    if (ErrorState())
        return 0;

    MarkInitialSeamsAsFaux(shell, baseMesh);

    int numPlacedCones = 0;

    for (;;) {
        std::vector<int> coneIndices;
        FindConesWithThreshold(cst, coneIndices);

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
                LOG_DEBUG << "Selected face " << tri::Index(shell, p.F()) << " edge " << p.E();
                tri::UpdateFlags<Mesh>::FaceClearFaceEdgeS(shell);
                PosF boundaryPos = SelectShortestSeamPathToBoundary(shell, p);
                //SelectShortestSeamPathToPeak(shell, p);
                RectifyCut(shell, boundaryPos);
                tri::CutMeshAlongSelectedFaceEdges(shell);
                CleanupShell(shell);
                tri::UpdateTopology<Mesh>::FaceFace(shell);
            }
            numPlacedCones += (int) coneIndices.size();
        } else {
            break;
        }
    }

    MarkInitialSeamsAsFaux(shell, baseMesh);
    if (numPlacedCones > 0) {
        // need to update the boundary information, which is implicitly computed
        // by the hole filling function (ugly...)
        ClearHoleFillingFaces(shell, true, true);
        CloseMeshHoles(shell);
    }

    return numPlacedCones;
}


/* What this method does is:
 *  1. Remove any hole filling faces from the shell
 *  2. Bring the shell in its 3D configuration
 *  3. Close the holes
 *  4. Iteratively place cones and cut until the maximum conformal scaling
 *     factor is below the 'cst' parameter. Cones are placed as in the paper
 *     from Ben-Chen et al. 'Conformal flattening by curvature prescription and
 *     metric scaling' (2008), restricted to the vertices that lie on uv seams
 *  6. Restore the UV configuration of the shell by reinitializing the solution
 *     and rebuilding the scaffold if needed */
int ParameterizerObject::PlaceCutWithConesUntilThreshold(double cst)
{
    if (ErrorState())
        return 0;

    // Put the shell in a configuration suitable to compute the conformal scaling factors
    ClearHoleFillingFaces(shell, true, true);

    SyncShellWith3D(shell);
    tri::UpdateTopology<Mesh>::FaceFace(shell);

    CloseMeshHoles(shell);

    int numPlacedCones = PlaceCutWithConesUntilThreshold_3D(cst);

    // if cones were found we reinitialize the solution, otherwise just leave
    // everything as is. In this second case we can simply return the shell to
    // the previous uv state, which is not changed.

    if (numPlacedCones > 0) {
        InitializeSolution();
    }

    // Cleanup (and restore hole-filling and scaffold faces if needed)
    ClearHoleFillingFaces(shell, true, true);
    SyncShellWithUV(shell);

    if (strategy.padBoundaries)
        CloseShellHoles(shell, strategy.geometry, baseMesh);

    if (strategy.scaffold) {
        BuildScaffold(shell, strategy.geometry, baseMesh);
    }

    MarkInitialSeamsAsFaux(shell, baseMesh);

    if (OptimizerIsInitialized()) {
        opt->UpdateCache();
    }

    return numPlacedCones;
}

void ParameterizerObject::FindCones(int ncones, std::vector<int>& coneIndices)
{
    coneIndices.clear();

    ComputeDistanceFromBorderOnSeams(shell);

    for (int i = 0; i < ncones; ++i) {
        Eigen::VectorXd csf;
        if (!ComputeConformalScalingFactors(csf, coneIndices)) {
            LOG_WARN << "Scaling factors computation failed";
            return;
        }

        // Select seam vertices
        tri::UpdateFlags<Mesh>::VertexClearS(shell);
        for (auto& sv : shell.vert)
            if (sv.Q() < Infinity())
                sv.SetS();

        // Find the max scale factor at seam vertices
        double maxscale = 0.0;
        double minscale = Infinity();
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

        LOG_DEBUG << "CONFORMAL SCALING FACTORS: max = " << maxscale << " | min = " << minscale <<  " | ratio = " << minscale / maxscale;

        coneIndices.push_back(max_k);
    }
}

/* This function is a bit messed up, in theory one could choose an arbitrary number
 * of cone vertices, but it makes little sense to me to treat a 'designated' cone
 * as boundary vertex without actually cutting the surface. So I restrict the number
 * of cones to 1, forcing to cut the surface repeatedly. */
void ParameterizerObject::FindConesWithThreshold(double cst, std::vector<int>& coneIndices)
{
    constexpr int MAX_NUM_CONES = 1;

    coneIndices.clear();

    ComputeDistanceFromBorderOnSeams(shell);
    bool seamsToBoundary = false;
    for (auto& sv : shell.vert) {
        if (!sv.IsB() && sv.Q() < Infinity()) {
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
            LOG_WARN << "Warning: scaling factors computation failed";
            return;
        }

        // Select seam vertices that reach the boundary
        tri::UpdateFlags<Mesh>::VertexClearS(shell);
        for (auto& sv : shell.vert)
            if (sv.Q() > 0 && sv.Q() < Infinity())
                sv.SetS();

        // Find the max scale factor at seam vertices
        double maxscale = 0.0;
        double minscale = Infinity();
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

        LOG_DEBUG << "CONFORMAL SCALING FACTORS: max = " << maxscale << " | min = " << minscale <<  " | ratio = " << minscale / maxscale;

        if (maxscale > cst)
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

    LOG_VERBOSE << "Computing conformal scaling factors...";

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

    // Compute the metric scaling of the new target curvature
    ComputeCotangentWeightedLaplacian(m, M, index);
    M += Eigen::VectorXd::Constant(M.rows(), 1e-12).asDiagonal();

    Eigen::SparseMatrix<double> M_internal = M.block(0, 0, n_internal, n_internal);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverLDLT;
    solverLDLT.compute(M_internal);
    if (solverLDLT.info() != Eigen::Success) {
        LOG_WARN << "ComputeConformalScalingFactors(): preliminary factorization failed";
        return false;
    }

    // Compute the new curvature vector Knew
    Eigen::VectorXd Knew = Eigen::VectorXd::Zero(m.VN());
    Eigen::VectorXd G(n_internal);

    // This loops takes a lot of time (a backward solve for each boundary vertex...)
    for (int i = 0; i < n_boundary; ++i) {
        G = solverLDLT.solve(-M.block(0, n_internal+i, n_internal, 1));
        if (solverLDLT.info() != Eigen::Success) {
            LOG_WARN << "ComputeConformalScalingFactors(): backward substitution to compute G(" << i << ") failed";
            return false;
        }
        Knew[n_internal+i] = Korig[n_internal+i] + G.dot(Korig.head(n_internal));
    }

    LOG_DEBUG << "Computation of the target curvatures took " << timer.TimeSinceLastCheck() << " seconds";

    solverLDLT.compute(M);
    if (solverLDLT.info() != Eigen::Success) {
        LOG_WARN << "ComputeConformalScalingFactors(): factorization of the cotangent-weighted Laplacian failed";
        return false;
    }

    Eigen::VectorXd p = Knew - Korig;

    Eigen::VectorXd scalingFactors = solverLDLT.solve(Knew - Korig);
    csf.resize(scalingFactors.size());
    for (int i = 0; i < m.VN(); ++i) {
        csf[i] = scalingFactors[index[m.vert[i]]];
    }

    double maxscale = 0.0;
    double minscale = Infinity();
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

    LOG_DEBUG << "Solution of the scaling factors system took " << timer.TimeSinceLastCheck() << " seconds";
    LOG_DEBUG << "Overall time: " << timer.TimeElapsed() << " seconds";

    return true;
}


