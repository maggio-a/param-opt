#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/export.h>

#include <string>
#include <vector>
#include <iostream>

#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>

#include "mesh.h"
#include "mesh_graph.h"
#include "uv.h"
#include "texture_optimization.h"
#include "texture_rendering.h"
#include "timer.h"
#include "gl_utils.h"

using namespace vcg;

void LogStrategy(ParameterizationStrategy strategy, double tol)
{
    std::cout << "[LOG] Parameterization strategy: ";
    std::cout << "iterations=" << strategy.optimizerIterations << " , "
              << "tolerance=" << tol << " "
              << "padded inner boundaries=" << strategy.padBoundaries
              << std::endl;
}

void LogDistortionStats(std::shared_ptr<MeshGraph> graph)
{
    vcg::Distribution<double> d;

    std::cout << "[LOG] Distortion" << std::endl;
    std::cout << "[LOG] Type Source Min Max Avg Variance" << std::endl;

    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Area, ParameterizationGeometry::Texture);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Area Texture " << d.Min() << " " << d.Max() << " " << d.Avg() << " " << d.Variance() << std::endl;


    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Angle, ParameterizationGeometry::Texture);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Angle Texture " << d.Min() << " " << d.Max() << " " << d.Avg() << " " << d.Variance() << std::endl;


    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Area, ParameterizationGeometry::Model);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Area Model " << d.Min() << " " << d.Max() << " " << d.Avg() << " " << d.Variance() << std::endl;


    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Angle, ParameterizationGeometry::Model);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Angle Model " << d.Min() << " " << d.Max() << " " << d.Avg() << " " << d.Variance() << std::endl;

}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " model [minRegionSize(int)] [--nofilter]" << std::endl;
        return -1;
    }

    bool filter = true;
    if (argc > 3 && std::string("--nofilter").compare(argv[3]) == 0) {
        std::cout << "Pull-Push filter disabled" << std::endl;
        filter = false;
    } else {
        std::cout << "Pull-Push filter enabled" << std::endl;
    }

    int minRegionSize = -1;
    if (argc > 2) minRegionSize = atoi(argv[2]);

    Mesh m;
    TextureObjectHandle textureObject;
    int loadMask;
    std::string modelName;

    if (LoadMesh(m, argv[1], textureObject, loadMask, modelName) == false) {
        std::cout << "Failed to open mesh" << std::endl;
        std::exit(-1);
    }

    tri::UpdateTopology<Mesh>::FaceFace(m);
    if (tri::Clean<Mesh>::CountNonManifoldEdgeFF(m, false)) {
        std::cout << "Mesh is not edge manifold" << std::endl;
        std::exit(-1);
    }

    if (minRegionSize == -1) {
        if (m.FN() < 100000) minRegionSize = 1000;
        else if (m.FN() < 300000) minRegionSize = 5000;
        else minRegionSize = 10000;
    }

    assert(m.FN() < 2000000);

    assert(textureObject->ArraySize() == 1 && "Currently only single texture is supported");

    assert(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    if (minRegionSize > m.FN()) {
        std::cout << "WARNING: minFaceCount > m.FN()" << std::endl;
    }


    ParameterizationStrategy strategy;
    strategy.directParameterizer = FixedBorderBijective;
    strategy.optimizer = SymmetricDirichletOpt;
    strategy.geometry = Texture;
    strategy.descent = ScalableLocallyInjectiveMappings;
    strategy.optimizerIterations = 600;
    strategy.padBoundaries = true;

    double tolerance = 0.0001;

    PreprocessMesh(m);
    StoreWedgeTexCoordAsAttribute(m);

    LogStrategy(strategy, tolerance);

    float uvMeshBorder;
    auto graph = ComputeParameterizationGraph(m, textureObject, &uvMeshBorder);

    double areaUvBefore = graph->AreaUV();

    GLInit();

    // Print original info
    PrintParameterizationInfo(graph);

    Timer t;

    GraphManager gm{graph};

    std::cout << "[LOG] Edge weight limit (minRegionSize) = " << minRegionSize << std::endl;
    ReduceTextureFragmentation_NoPacking(gm, minRegionSize);

    int c = ParameterizeGraph(gm, strategy, true, tolerance, true);
    if (c > 0) std::cout << "WARNING: " << c << " regions were not parameterized correctly" << std::endl;

    std::cout << "[LOG] Scale factor of the packed chart = " << graph->AreaUV() / areaUvBefore << std::endl;

    LogDistortionStats(graph);

    std::cout << "Rendering texture..." << std::endl;
    TextureObjectHandle newTexture = RenderTexture(m, textureObject, filter, nullptr);

    std::cout << "Processing took " << t.TimeElapsed() << " seconds" << std::endl;

    if (SaveMesh(m, modelName.c_str(), newTexture) == false) {
        std::cout << "Model not saved correctly" << std::endl;
    }

    auto graph2 = ComputeParameterizationGraph(m, textureObject, &uvMeshBorder);
    // Print optimized info
    PrintParameterizationInfo(graph2);

    GLTerminate();

    return 0;
}


