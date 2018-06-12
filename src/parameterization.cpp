#include "parameterization.h"
#include "mesh_graph.h"
#include "energy.h"
#include "iterative.h"
#include "parameterization_checker.h"
#include "mesh_utils.h"
#include "timer.h"
#include "uniform_solver.h"

#include <memory>

#include <vcg/complex/algorithms/crease_cut.h>

#include <Eigen/Core>

#include <wrap/io_trimesh/export.h>


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
      energyDiffTolerance{1e-6}
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
    assert(init && "Failed to initialize ParameterizerObject solution");
    assert(CheckUVConnectivity(shell));
    assert(CheckLocalInjectivity(shell));

    if (strategy.padBoundaries == false)
        ClearHoleFillingFaces(shell);
    InitializeOptimizer();
    SyncShellWithUV(shell);
    MarkInitialSeamsAsFaux(shell, baseMesh);
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
    std::vector<int> coneIndices = {};
    assert(ComputeConformalScalingFactors(csf, coneIndices));
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


bool ParameterizerObject::InitializeSolution()
{
    UniformSolver<Mesh> u{shell};
    bool solved = u.Solve();
    if (solved) SyncShellWithUV(shell);
    return solved;
}

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
    SyncShellWithUV(shell);
    iterationCount++;
    return info;
}

void ParameterizerObject::RemeshHolefillingAreas()
{
    RemeshShellHoles(shell, strategy.geometry, baseMesh);
    MarkInitialSeamsAsFaux(shell, baseMesh);
    opt->UpdateCache();
}

void ParameterizerObject::PlaceCut()
{
    ComputeDistanceFromBorderOnSeams(shell);

    double maxDistance = 0;
    for (auto& sv : shell.vert) {
        if (sv.Q() < 0) assert(0);
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

                double weightedEnergy = energy->E(sf, false) * w;
                // w is INFINITY if the seam does not reach the mesh boundary
                if (std::isfinite(weightedEnergy) && weightedEnergy > maxEnergy) {
                    maxEnergy = weightedEnergy;
                    startFace = &sf;
                    cutEdge = i;
                }
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

/* What this method does is:
 *  1. Remove any hole filling faces from the shell
 *  2. Sync the shell with the model
 *  3. Close the holes
 *  4. Compute 'ncones' cone singularities using the method described in
 *     Ben-Chen et al. 'Conformal flattening by curvature prescription and
 *     metric scaling' (2008), restricted to the vertices that lie on uv seams
 *  5. Iteratively cut the shell, connecting each cone with the boundary
 *  6. Remove the hole filling faces, sync the shell with its UV configuration
 *     and close the holes in UV space (this is necessary to avoid potential
 *     flips) */
bool ParameterizerObject::PlaceCutWithConeSingularities(int ncones)
{
    assert(ncones > 0);

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
    ClearHoleFillingFaces(shell);
    SyncShellWithModel(shell, baseMesh);
    for (auto& sf : shell.face) {
        for (int i = 0; i < 3; ++i) {
            if (sf.IsF(i)) {
                sf.C() = vcg::Color4b::Blue;
                sf.FFp(i)->C() = vcg::Color4b::Blue;
            }
        }
    }
    //tri::io::Exporter<Mesh>::Save(shell, "shell_clear.obj", tri::io::Mask::IOM_ALL);
    ComputeBoundaryInfo(shell);
    CloseMeshHoles(shell);
    MarkInitialSeamsAsFaux(shell, baseMesh);

    std::vector<int> coneIndices;
    FindCones(ncones, coneIndices, true);

    ComputeDistanceFromBorderOnSeams(shell);

    // At this point the cones are guaranteed to be placed on a mesh vertex
    tri::UpdateFlags<Mesh>::VertexClearS(shell);
    for (int i : coneIndices)
        shell.vert[i].SetS();

    std::vector<PosF> pv;
    for (auto& sf : shell.face) {
        if (!sf.holeFilling) {
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
    }

    ClearHoleFillingFaces(shell);
    ComputeBoundaryInfo(shell); // boundary changed after the cut
    SyncShellWithUV(shell);
    CloseShellHoles(shell, strategy.geometry, baseMesh);
    RemeshShellHoles(shell, strategy.geometry, baseMesh);

    if (OptimizerIsInitialized())
        opt->UpdateCache();

    //tri::io::Exporter<Mesh>::Save(shell, "shell_after_cut.obj", tri::io::Mask::IOM_ALL);
    return true;
}

void ParameterizerObject::FindCones(int ncones, std::vector<int>& coneIndices, bool onSeams)
{
    coneIndices.clear();

    if (onSeams)
        ComputeDistanceFromBorderOnSeams(shell);

    for (int i = 0; i < ncones; ++i) {
        Eigen::VectorXd csf;
        ComputeConformalScalingFactors(csf, coneIndices);

        // Select seam vertices
        if (onSeams) {
            tri::UpdateFlags<Mesh>::VertexClearS(shell);
            for (auto& sv : shell.vert)
                if (sv.Q() < INFINITY)
                    sv.SetS();
        } else {
            tri::UpdateFlags<Mesh>::VertexSetS(shell);
        }

        // Find the max scale factor as seam vertices
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

        coneIndices.push_back(max_k);
    }
}

void ParameterizerObject::ProbeCut()
{

}

using Td = Eigen::Triplet<double>;

/* assumes that cone singularity vertices have the selected flag set */
void ComputeMarkovMatrix(Mesh& m, Eigen::SparseMatrix<double>& L, SimpleTempData<Mesh::VertContainer, int>& Idx)
{
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

void ComputeCotangentWeightedLaplacian(Mesh& shell, Mesh& baseMesh, Eigen::SparseMatrix<double>& L,
                                       SimpleTempData<Mesh::VertContainer, int>& Idx,
                                       Mesh::PerFaceAttributeHandle<int>& ia)
{
    L.resize(shell.VN(), shell.VN());
    L.setZero();
    std::vector<Td> tri;
    double eps = std::numeric_limits<double>::epsilon();
    for (auto &sf : shell.face) {
        for (int i = 0; i < 3; ++i) {
            Mesh::VertexPointer vi = sf.V0(i);
            Mesh::VertexPointer vj = sf.V1(i);
            Mesh::VertexPointer vk = sf.V2(i);

            assert(Idx[vi] >= 0 && Idx[vi] < shell.VN());
            assert(Idx[vj] >= 0 && Idx[vj] < shell.VN());
            assert(Idx[vk] >= 0 && Idx[vk] < shell.VN());

            double alpha_ij = std::max(VecAngle(vi->P() - vk->P(), vj->P() - vk->P()), eps); // angle at k
            double alpha_ik = std::max(VecAngle(vi->P() - vj->P(), vk->P() - vj->P()), eps); // angle at j
            double weight_ij = 0.5 * std::tan(M_PI_2 - alpha_ij);
            double weight_ik = 0.5 * std::tan(M_PI_2 - alpha_ik);
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
    //ClearHoleFillingFaces(shell);

    /*
    Mesh m;
    CopyToMesh(*chart, m);
    tri::UpdateTopology<Mesh>::FaceFace(m);
    ComputeBoundaryInfo(m);
    CloseMeshHoles(m);

    tri::Clean<Mesh>::SplitNonManifoldVertex(m, 0.0);
    */
    Mesh& m = shell;

    //vcg::tri::io::Exporter<Mesh>::Save(m, "shell_model.ply");

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


    /*
    ComputeMarkovMatrix(m, M, index);

    Eigen::SparseMatrix<double> Q{n_internal, n_internal};
    Q.setIdentity();
    Q -= M.block(0, 0, n_internal, n_internal);

    Eigen::SparseMatrix<double> T = M.block(0, n_internal, n_internal, n_boundary);

    Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>, Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::StorageIndex>> solver;
    solver.compute(Q);
    if (solver.info() != Eigen::Success) {
        std::cout << "Factorization of the Markov sub-matrix failed" << std::endl;
        return false;
    }

    std::cout << "Factorization of the Markov matrix took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    */

    // Compute the new curvature vector Knew
    // Note that this is computed in the original model space
    auto ia = GetFaceIndexAttribute(m);

    Eigen::VectorXd Korig{m.VN()};
    Korig.setZero();
    for (auto& sf : m.face) {
        //auto& f = baseMesh.face[ia[sf]];
        for (int i = 0; i < 3; ++i) {
            Korig[index[sf.V(i)]] -= VecAngle(sf.P1(i) - sf.P0(i), sf.P2(i) - sf.P0(i));
        }
    }
    Korig.segment(0, n_internal) += Eigen::VectorXd::Constant(n_internal, 2.0 * M_PI);
    Korig.segment(n_internal, n_boundary) += Eigen::VectorXd::Constant(n_boundary, M_PI);

    std::cout << "Sum Korig = " << Korig.sum() << std::endl;

    Eigen::VectorXd Knew{m.VN()};
    Knew.setZero();
    Eigen::VectorXd G{n_internal};





    /*
    for (int i = 0; i < n_boundary; ++i) {
        G = solver.solve(T.col(i));
        if (solver.info() != Eigen::Success) {
            std::cout << "Backward substitution to compute G(" << i << ") failed" << std::endl;
            return false;
        }
        Knew[n_internal+i] = Korig[n_internal+i] + G.dot(Korig.head(n_internal));
    }
    std::cout << "Sum Knew = " << Knew.sum() << std::endl;
    std::cout << "Computation of the target curvatures took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    */






    // Compute the metric scaling of the new target curvature
    ComputeCotangentWeightedLaplacian(m, baseMesh, M, index, ia);
    //M += Eigen::VectorXd::Constant(M.rows(), 1e-12).asDiagonal();




    Eigen::SparseMatrix<double> L = M.block(0, 0, n_internal, n_internal);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> CHOL;
    CHOL.compute(L);
    if (CHOL.info() != Eigen::Success) {
        assert(0);
    }
    Eigen::VectorXd delta{n_internal};
    delta.setZero();
    for (int i = 0; i < n_boundary; ++i) {
        delta[i] = 1.0;
        //G = CHOL.solve(T.col(i));
        G = CHOL.solve(-M.block(0, n_internal+i, n_internal, 1));
        //G = CHOL.solve(delta);
        if (CHOL.info() != Eigen::Success) {
            std::cout << "Backward substitution to compute G(" << i << ") failed" << std::endl;
            return false;
        }
        delta[i] = 0.0;
        Knew[n_internal+i] = Korig[n_internal+i] + G.dot(Korig.head(n_internal));
    }
    std::cout << "Sum Knew = " << Knew.sum() << std::endl;
    std::cout << "Computation of the target curvatures took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;



    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholesky;
    cholesky.compute(M);
    if (cholesky.info() != Eigen::Success) {
        std::cout << "Factorization of the cotangent-weighted Laplacian failed" << std::endl;
        return false;
    }

    std::cout << "Factorization of the cotangent-weighted Laplacian took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;

    Eigen::VectorXd scalingFactors = cholesky.solve(Knew - Korig);
    csf.resize(scalingFactors.size());
    for (int i = 0; i < m.VN(); ++i) {
        csf[i] = scalingFactors[index[m.vert[i]]];
    }

    double maxscale = 0.0;
    double minscale = std::numeric_limits<double>::infinity();
    int max_i = -1;
    int min_i = -1;
    for (int i = 0; i < m.VN(); ++i) {
        double s = std::abs(csf[i]);
        //assert(scalingFactors[i] > 0.0);
        if (s > maxscale) {
            maxscale = s;
            max_i = i;
        }
        if (s < minscale) {
            minscale = s;
            min_i = i;
        }
    }

    std::cout << " max scale at vertex " << max_i << std::endl;
    std::cout << " min scale at vertex " << min_i << std::endl;

    //coneIdx.push_back(max_i);

    std::cout << "Solution of the scaling factors system took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    std::cout << "Overall time: " << timer.TimeElapsed() << " seconds" << std::endl;

    return true;
}







