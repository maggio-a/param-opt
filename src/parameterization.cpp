#include "parameterization.h"
#include "mesh_graph.h"
#include "energy.h"
#include "iterative.h"
#include "parameterization_checker.h"
#include "mesh_utils.h"
#include "timer.h"

#include <memory>

#include <vcg/complex/algorithms/crease_cut.h>

#include <Eigen/Core>



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
      energyDiffTolerance{1e-9},
      stats{0}
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

    /*
     * In case of warm start, we should build the shell according to the existing
     * parameterization (after having checked that is valid wrt the given strategy
     * */
    std::cout << "TODO FIXME WARM START CODE" << std::endl;
    BuildShell(shell, *chart, strategy.geometry);
    std::cout << "WARNING forcing warm start to true" << std::endl;
    strategy.warmStart = true;
    if (strategy.padBoundaries == false)
        ClearHoleFillingFaces(shell);

    InitializeOptimizer();
    SyncShell(shell);

    MarkInitialSeamsAsFaux(shell, baseMesh);
}

void ParameterizerObject::Sync()
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
    Eigen::MatrixXd D = opt->ComputeDescentDirection();
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
}

void ParameterizerObject::MapConformalScalingFactorsToShellVertexColor()
{
    Eigen::VectorXd csf;
    assert(ComputeConformalScalingFactors(csf));
    assert(csf.size() == shell.VN());
    for (int i = 0; i < shell.VN(); ++i) {
        shell.vert[i].Q() = csf[i];
    }
    tri::UpdateColor<Mesh>::PerFaceQualityRamp(shell);
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
    default:
        assert(0);
    }
    energy = energy_sd;
}

bool ParameterizerObject::OptimizerIsInitialized()
{
    return energy && opt;
}


/* Parameterization-related methods
 * ================================ */


bool ParameterizerObject::Parameterize()
{
    if (!OptimizerIsInitialized())
        InitializeOptimizer();

    if (strategy.warmStart) {
        assert(CheckLocalInjectivity(shell));
        assert(CheckUVConnectivity(shell));
    }

    // check if remeshing is required during iterations
    needsRemeshing = false;
    for (auto& sf : shell.face) {
        if (sf.holeFilling) needsRemeshing = true;
        double areaUV = (sf.V(1)->T().P() - sf.V(0)->T().P()) ^ (sf.V(2)->T().P() - sf.V(0)->T().P());
        assert(areaUV > 0 && "Parameterization is not bijective");
    }

    // parameterizer state
    Timer t;
    int i;

    IterationInfo info;
    for (i = 0; i < strategy.optimizerIterations; ++i) {
        info = Iterate();
        if (info.gradientNorm < gradientNormTolerance) {
            std::cout << "Stopping because gradient magnitude is small enough (" << info.gradientNorm << ")" << std::endl;
            break;
        }
        if (info.energyDiff < energyDiffTolerance) {
            std::cout << "Stopping because energy improvement is too small (" << info.energyDiff << ")" << std::endl;
            break;
        }
        if (needsRemeshing && (i > 0) && (i % 30) == 0)
            RemeshHolefillingAreas();
    }
    std::cout << "Stopped after " << i << " iterations, gradient magnitude = " << info.gradientNorm
              << ", energy value = " << energy->E_IgnoreMarkedFaces() << std::endl;

    float minq, maxq;
    tri::Stat<Mesh>::ComputePerFaceQualityMinMax(shell, minq, maxq);
    std::cout << "Min distortion value = " << minq << std::endl;
    std::cout << "Max distortion value = " << maxq << std::endl;

    std::cout << "Optimization took " << t.TimeSinceLastCheck() << " seconds" << std::endl;

    return true;
}

IterationInfo ParameterizerObject::Iterate()
{
    if (!OptimizerIsInitialized())
        InitializeOptimizer();
    IterationInfo info;
    info.energyVal = opt->Iterate(info.gradientNorm, info.energyDiff);
    SyncShell(shell);
    iterationCount++;
    return info;
}

void ParameterizerObject::RemeshHolefillingAreas()
{
    RemeshShellHoles(shell, strategy.geometry, baseMesh);
    MarkInitialSeamsAsFaux(shell, baseMesh);
    opt->UpdateCache();
}

//#include <wrap/io_trimesh/export.h>
void ParameterizerObject::PlaceCut()
{
    double maxDistance = ComputeDistanceFromBorderOnSeams(shell);

    MapEnergyToShellFaceColor();

    // Select candidate crease edge for cutting
    double maxEnergy = 0;
    Mesh::FacePointer startFace = nullptr;
    int cutEdge = -1;
    for (auto& sf : shell.face) {
        for (int i = 0; i < 3; ++i) if (sf.IsF(i)) {
            double w = std::max(sf.V0(i)->Q(), sf.V1(i)->Q()) / maxDistance;

            if (w < INFINITY) w = 1.0;

            double weightedEnergy = energy->E(sf, false) * w;
            // w is INFINITY if the seam does not reach the mesh boundary
            if (std::isfinite(weightedEnergy) && weightedEnergy > maxEnergy) {
                maxEnergy = weightedEnergy;
                startFace = &sf;
                cutEdge = i;
            }
        }
    }
    assert(startFace != nullptr);

    //tri::io::Exporter<Mesh>::Save(shell, "shell_before_cut.obj", tri::io::Mask::IOM_ALL);

    std::cout << "Selected face " << tri::Index(shell, startFace) << " edge " << cutEdge << std::endl;

    PosF p(startFace, cutEdge);
    assert(!p.IsBorder());

    tri::UpdateFlags<Mesh>::FaceClearFaceEdgeS(shell);
    PosF boundaryPos = SelectShortestSeamPathToBoundary(shell, p);
    //SelectShortestSeamPathToPeak(shell, p);
    bool corrected = RectifyCut(shell, boundaryPos);

    tri::CutMeshAlongSelectedFaceEdges(shell);
    CleanupShell(shell);

    tri::UpdateTopology<Mesh>::FaceFace(shell);

    opt->UpdateCache();
    //tri::io::Exporter<Mesh>::Save(shell, "shell_after_cut.obj", tri::io::Mask::IOM_ALL);
}

void ParameterizerObject::ProbeCut()
{

}

using Td = Eigen::Triplet<double>;

/* assumes that cone singularity vertices have the selected flag set */
void ComputeMarkovMatrix(Mesh& m, Eigen::SparseMatrix<double>& L, SimpleTempData<Mesh::VertContainer, int>& Idx)
{
    L.resize(m.VN(), m.VN());
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
            for (int i = 0; i < 3; ++i) {
                double c = entry.second.size();
                for (auto & vp : entry.second) {
                    tri.push_back(Td(Idx[entry.first], Idx[vp], (1.0 / c)));
                }
            }
        } else {
            tri.push_back(Td(Idx[entry.first], Idx[entry.first], 1.0));
        }
    }
    L.setFromTriplets(tri.begin(), tri.end());
    L.makeCompressed();
}

void ComputeCotangentWeightedLaplacian(Mesh& shell, Mesh& baseMesh, Eigen::SparseMatrix<double>& L,
                                       SimpleTempData<Mesh::VertContainer, int>& Idx,
                                       Mesh::PerFaceAttributeHandle<int>& ia)
{
    L.resize(shell.VN(), shell.VN());
    std::vector<Td> tri;
    double eps = std::numeric_limits<double>::epsilon();
    for (auto &sf : shell.face) {
        auto& f = baseMesh.face[ia[sf]];
        for (int i = 0; i < 3; ++i) {
            Mesh::VertexPointer vi = f.V0(i);
            Mesh::VertexPointer vj = f.V1(i);
            Mesh::VertexPointer vk = f.V2(i);
            double alpha_ij = std::max(VecAngle(vi->P() - vk->P(), vj->P() - vk->P()), eps); // angle at k
            double alpha_ik = std::max(VecAngle(vi->P() - vj->P(), vk->P() - vj->P()), eps); // angle at j
            double weight_ij = std::tan(M_PI_2 - alpha_ij);
            double weight_ik = std::tan(M_PI_2 - alpha_ik);
            tri.push_back(Td(Idx[vi], Idx[vj], weight_ij));
            tri.push_back(Td(Idx[vi], Idx[vk], weight_ik));
            tri.push_back(Td(Idx[vi], Idx[vi], -(weight_ij + weight_ik)));
        }
    }
    L.setFromTriplets(tri.begin(), tri.end());
    L.makeCompressed();
}

/* Reference: BenChen 2009 curvature prescription metric scaling */
bool ParameterizerObject::ComputeConformalScalingFactors(Eigen::VectorXd& csf)
{
    ClearHoleFillingFaces(shell);

    Timer timer;

    std::cout << "Computing conformal scaling factors..." << std::endl;

    // Permutate indices so that singular vertices are the last ones
    SimpleTempData<Mesh::VertContainer, int> index{shell.vert};
    tri::UpdateFlags<Mesh>::VertexBorderFromFaceAdj(shell);
    tri::UpdateFlags<Mesh>::VertexClearS(shell);
    int n_internal = 0;
    int j = shell.VN() - 1;
    for (auto& v : shell.vert) {
        if (v.IsB()) {
            index[v] = j--;
            v.SetS();
        } else
            index[v] = n_internal++;
    }
    int n_boundary = shell.VN() - n_internal;

    // Compute Laplacian and Markov system
    Eigen::SparseMatrix<double> M;
    ComputeMarkovMatrix(shell, M, index);

    Eigen::SparseMatrix<double> Q{n_internal, n_internal};
    Q.setIdentity();
    Q -= M.block(0, 0, n_internal, n_internal);

    Eigen::SparseMatrix<double> T = M.block(0, n_internal, n_internal, n_boundary);

    Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>, Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::StorageIndex>> solver;
    solver.compute(T);
    if (solver.info() != Eigen::Success) {
        std::cout << "Factorization of the Markov sub-matrix failed" << std::endl;
        return false;
    }

    std::cout << "Factorization of the Markov matrix took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;

    Eigen::SparseMatrix<double> G{n_internal, n_internal};
    Eigen::VectorXd t{n_internal};
    t.setZero();
    for (int i = 0; i < n_boundary; ++i) {
        t[i] = 1.0;
        G.col(i) = solver.solve(t);
        if (solver.info() != Eigen::Success) {
            std::cout << "Backward substitution to compute G(" << i << ") failed" << std::endl;
            return false;
        }
        t[i] = 0.0;
    }

    std::cout << "Computation of G took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;

    // Compute the new curvature vector Knew
    // Note that this is computed in the original model space
    auto ia = GetFaceIndexAttribute(shell);

    Eigen::VectorXd Korig{shell.VN()};
    Eigen::VectorXd Knew{shell.VN()};
    Korig.setZero();
    Knew.setZero();

    for (auto& sf : shell.face) {
        auto& f = baseMesh.face[ia[sf]];
        for (int i = 0; i < 3; ++i) {
            Korig[index[sf.V(i)]] += VecAngle(f.P1(i) - f.P0(i), f.P2(i) - f.P0(i));
        }
    }
    Knew.tail(n_boundary) = Korig.tail(n_boundary) + G.transpose() * Korig.head(n_internal);

    std::cout << "Computation of the target curvatures took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;

    // Compute the metric scaling of the new target curvature
    ComputeCotangentWeightedLaplacian(shell, baseMesh, M, index, ia);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholesky;
    cholesky.compute(M);
    if (cholesky.info() != Eigen::Success) {
        std::cout << "Factorization of the cotangent-weighted Laplacian failed" << std::endl;
        return false;
    }

    std::cout << "Factorization of the cotangent-weighted Laplacian took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;

    Eigen::VectorXd scalingFactors = cholesky.solve(Korig - Knew);
    csf.resize(scalingFactors.size());
    for (int i = 0; i < shell.VN(); ++i) {
        csf[i] = scalingFactors[index[shell.vert[i]]];
    }

    std::cout << "Solution of the scaling factors system took" << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    std::cout << "Overall time: " << timer.TimeElapsed() << " seconds" << std::endl;

    return true;
}







