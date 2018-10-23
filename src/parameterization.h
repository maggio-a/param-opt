#ifndef PARAMETERIZATION_H
#define PARAMETERIZATION_H

#include "mesh.h"
#include "mesh_graph.h"
#include "energy.h"
#include "iterative.h"

#include <Eigen/Core>
#include <memory>


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


/* Builds a shell for the given chart. A shell is a mesh object specifically
 * constructed to compute the parameterization of a chart. In order to support
 * the various operations that we need to perform on it, it is more convenient
 * to keep its shape as the current 2D parameter-space configuration (possibly
 * not updated). The shell has suitable attributes to retrieve information about
 * the shell-face to input mesh-face mappings, as well as the target shape
 * features of each face to guide the parameterization process. (See also the
 * comments in mesh_attribute.h).
 * The shell is initialized using Tutte's parameterization, and it is scaled so
 * that it matches the target area. */
bool BuildShell(Mesh& shell, FaceGroup& fg, ParameterizationGeometry targetGeometry, bool useExistingUV);

/* This function synchronizes a shell with its UV coordinates, that is it
 * updates its vertex coordinates to match the parameter space configurations
 * (with z = 0). The operation is performed per-vertex. */
void SyncShellWithUV(Mesh& shell);

/* This function synchronizes a shell with the model space coordinates of its
 * chart. */
void SyncShellWith3D(Mesh& shell);

/* Removes any hole-filling face from the shell and compacts its containers */
void ClearHoleFillingFaces(Mesh& shell, bool holefill, bool scaffold);

/* Closes shell holes, updating the shell attributes accordingly. Note that
 * faces that have a direct correspondence with the input mesh faces are not
 * touched by this procedure.
 * WARNING FIXME XXX as of now, this operation performed on the shell is not
 * guaranteed to always elect the same boundary as the one to be preserved
 * because it choses the longest one, but the boundary of the shell can change
 * during optimization */
void CloseShellHoles(Mesh& shell, ParameterizationGeometry targetGeometry, Mesh &inputMesh);

void BuildScaffold(Mesh& shell, ParameterizationGeometry targetGeometry, Mesh &inputMesh);

void RebuildScaffold(Mesh& shell, ParameterizationGeometry targetGeometry, Mesh &inputMesh);

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

public:

    enum Status { NoInit, Initialized, OK, Error_UnfeasibleShell, Error_InitFailed, Error_IterationFailed };

private:

    Mesh shell;
    Mesh& baseMesh;
    ChartHandle chart;
    ParameterizationStrategy strategy;
    std::shared_ptr<Energy> energy;
    std::shared_ptr<DescentMethod> opt;

    int iterationCount;

    double gradientNormTolerance;
    double energyDiffTolerance;
    double conformalScalingThreshold;

    double targetArea;

    Status status;

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

    Mesh& Shell();
    ChartHandle GetChart();
    int IterationCount();

    void MapEnergyToShellFaceColor();
    void MapEnergyGradientToShellVertexColor();
    void MapDescentDirectionToShellVertexColor();
    void MapLocalGradientVarianceToShellVertexColor();
    void MapConformalScalingFactorsToShellVertexColor();
    void ClearShellFaceColor();

    void SetGradientNormTolerance(double tol);
    void SetEnergyDiffTolerance(double tol);


    double GetGradientNormTolerance();
    double GetEnergyDiffTolerance();

    void Initialize();

    Status GetStatus();
    bool ErrorState();
    void PrintStatus();

private:

    bool InitializeSolution();

    void SetStatus(Status s);

    void InitializeOptimizer();
    bool OptimizerIsInitialized();

    int PlaceCutWithConesUntilThreshold_3D(double conformalScalingThreshold);

    bool ComputeConformalScalingFactors(Eigen::VectorXd& csf, const std::vector<int>& coneIndices);
    void FindCones(int ncones, std::vector<int>& coneIndices);
    void FindConesWithThreshold(double conformalScalingThreshold, std::vector<int>& coneIndices);
};

#endif // PARAMETERIZATION_H

