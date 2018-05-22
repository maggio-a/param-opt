#ifndef PARAMETERIZATION_H
#define PARAMETERIZATION_H

#include "mesh.h"
#include "mesh_graph.h"
#include "energy.h"
#include "iterative.h"

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

    ParameterizationStrategy() = default;
    ParameterizationStrategy(const ParameterizationStrategy&) = default;
};

inline ParameterizationStrategy DefaultStrategy()
{
    return ParameterizationStrategy{
        DirectParameterizer::FixedBorderBijective,
        EnergyType::SymmetricDirichlet,
        ParameterizationGeometry::Model,
        DescentType::ScalableLocallyInjectiveMappings,
        0, false, false, false
    };
}

inline ParameterizationStrategy MakeStrategy(DirectParameterizer directParameterizer,
                                             EnergyType energy,
                                             ParameterizationGeometry geometry,
                                             DescentType descent,
                                             int optimizerIterations,
                                             bool padBoundaries,
                                             bool applyCut,
                                             bool warmStart)
{
    return ParameterizationStrategy{
        directParameterizer,
        energy,
        geometry,
        descent,
        optimizerIterations,
        padBoundaries,
        applyCut,
        warmStart
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

public:

    ParameterizerObject(ChartHandle c, ParameterizationStrategy strat);
    ~ParameterizerObject();

    bool Parameterize();
    IterationInfo Iterate();
    void PlaceCut();
    void RemeshHolefillingAreas();
    void Sync();
    void Reset();

    Mesh& Shell();
    int IterationCount();

    void MapEnergyToShellFaceColor();
    void MapEnergyGradientToShellVertexColor();
    void MapDescentDirectionToShellVertexColor();
    void ClearShellFaceColor();

private:

    void InitializeOptimizer();
    bool OptimizerIsInitialized();
    void ProbeCut();
};

#endif // PARAMETERIZATION_H

