#ifndef PARAMETERIZATION_H
#define PARAMETERIZATION_H

#include "mesh.h"
#include "mesh_graph.h"
#include "energy.h"
#include "iterative.h"

#include <Eigen/Core>
#include <memory>

// TODO move to iterative.h
struct IterationInfo {
    double energyVal;
    double energyDiff;
    double gradientNorm;
};

enum DirectParameterizer {
    DCP, FixedBorderBijective, None
};

struct ParameterizationStrategy {
    DirectParameterizer directParameterizer;
    EnergyType energy;
    ParameterizationGeometry geometry;
    DescentType descent;
    int optimizerIterations;
    bool padBoundaries; // keep holes filled while optimizing from FixedBorderBijective
    bool applyCut;
    bool warmStart;
    bool scaffold;

    ParameterizationStrategy() = delete;
    ParameterizationStrategy(const ParameterizationStrategy&) = default;
};

inline ParameterizationStrategy DefaultStrategy()
{
    return ParameterizationStrategy{
        DirectParameterizer::FixedBorderBijective,
        EnergyType::SymmetricDirichlet,
        ParameterizationGeometry::Model,
        DescentType::ScalableLocallyInjectiveMappings,
        0, false, false, false, false
    };
}

inline ParameterizationStrategy MakeStrategy(DirectParameterizer directParameterizer,
                                             EnergyType energy,
                                             ParameterizationGeometry geometry,
                                             DescentType descent,
                                             int optimizerIterations,
                                             bool padBoundaries,
                                             bool applyCut,
                                             bool warmStart,
                                             bool scaffold)
{
    return ParameterizationStrategy{
        directParameterizer,
        energy,
        geometry,
        descent,
        optimizerIterations,
        padBoundaries,
        applyCut,
        warmStart,
        scaffold
    };
}

class ParameterizerObject {

    Mesh shell;
    Mesh& baseMesh;
    ChartHandle chart;
    ParameterizationStrategy strategy;
    std::shared_ptr<Energy> energy;
    std::shared_ptr<DescentMethod> opt;
    bool needsRemeshing;
    int iterationCount;

    double gradientNormTolerance;
    double energyDiffTolerance;

    double targetArea;

public:

    ParameterizerObject(ChartHandle c, ParameterizationStrategy strat);
    ~ParameterizerObject();

    bool Parameterize();
    IterationInfo Iterate();
    void PlaceCut();
    bool PlaceCutWithConeSingularities(int ncones);
    int PlaceCutWithConesUntilThreshold(double conformalScalingThreshold);
    void RemeshHolefillingAreas();

    /* Transfers the UV coordinates from the shell to the chart */
    void SyncChart();

    void Reset();
    bool InitializeSolution();

    Mesh& Shell();
    ChartHandle GetChart();
    int IterationCount();

    void MapEnergyToShellFaceColor();
    void MapEnergyGradientToShellVertexColor();
    void MapDescentDirectionToShellVertexColor();
    void MapLocalGradientVarianceToShellVertexColor();
    void MapConformalScalingFactorsToShellVertexColor();
    void ClearShellFaceColor();

    void ForceWarmStart();

    void SetGradientNormTolerance(double tol);
    void SetEnergyDiffTolerance(double tol);

    double GetGradientNormTolerance();
    double GetEnergyDiffTolerance();

private:

    void InitializeOptimizer();
    bool OptimizerIsInitialized();

    bool ComputeConformalScalingFactors(Eigen::VectorXd& csf, const std::vector<int>& coneIndices);
    void FindCones(int ncones, std::vector<int>& coneIndices);
    void FindConesWithThreshold(double conformalScalingThreshold, std::vector<int>& coneIndices);
};

#endif // PARAMETERIZATION_H

