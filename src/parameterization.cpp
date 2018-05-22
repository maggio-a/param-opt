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

ParameterizerObject::ParameterizerObject(ChartHandle c, ParameterizationStrategy strat)
    : shell{},
      baseMesh{c->mesh},
      chart{c},
      strategy(strat),
      energy{},
      opt{},
      needsRemeshing{false},
      iterationCount{0}
{
    Reset();
}

ParameterizerObject::~ParameterizerObject()
{
    // empty
}

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
        if (info.gradientNorm < 1e-3) {
            std::cout << "Stopping because gradient magnitude is small enough (" << info.gradientNorm << ")" << std::endl;
            break;
        }
        if (info.energyDiff < 1e-9) {
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

void ParameterizerObject::PlaceCut()
{
    double maxDistance = ComputeDistanceFromBorderOnSeams(shell);

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

            if (w < INFINITY) w = 1.0;

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
    //SelectShortestSeamPathToPeak(shell, p);

    tri::CutMeshAlongSelectedFaceEdges(shell);
    CleanupShell(shell);

    opt->UpdateCache();
}

void ParameterizerObject::MapEnergyToShellFaceColor()
{
    energy->MapToFaceQuality(false);
    tri::UpdateColor<Mesh>::PerFaceQualityRamp(shell, 0, 10.0);
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

void ParameterizerObject::ClearShellFaceColor()
{
    for (auto& sf : shell.face) {
        sf.C() = vcg::Color4b::White;
    }
}

void ParameterizerObject::Sync()
{
    for (std::size_t i = 0; i < chart->fpVec.size(); ++i) {
        for (int k = 0; k < 3; ++k) {
            chart->fpVec[i]->WT(k).P() = shell.face[i].WT(k).P();
        }
    }
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

void ParameterizerObject::ProbeCut()
{

}
