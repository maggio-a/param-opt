#include "mesh.h"
#include "mesh_graph.h"
#include "uv.h"
#include "texture_optimization.h"
#include "texture_rendering.h"
#include "timer.h"
#include "gl_utils.h"
#include "texture.h"
#include "mesh_viewer.h"
#include "logging.h"
#include "utils.h"

#include <wrap/io_trimesh/io_mask.h>

#include <string>
#include <vector>
#include <iostream>

#include <QCoreApplication>
#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>


using namespace vcg;


int MainCmd(Mesh& m, GraphHandle graph, TextureObjectHandle textureObject, Args args)
{
    ParameterizationStrategy strategy = MakeStrategy(
            DirectParameterizer::FixedBorderBijective,
            EnergyType::SymmetricDirichlet,
            ParameterizationGeometry::Texture,
            DescentType::CompositeMajorization,
            500,            // Number of iterations
            true,           // Fill holes
            true,           // Use cuts
            false,          // No warm start
            true            // Enable scaffold
    );
    double tolerance = 0.0002;

    LogExecutionParameters(args, strategy);

    /*
    LogStrategy(strategy, tolerance);

    std::vector<std::pair<int, Mesh::FacePointer>> cc;
    tri::Clean<Mesh>::ConnectedComponents(m, cc);
    int smallcomponents = 0;
    for (auto p : cc) if (p.first < 100) smallcomponents++;

    std::cout << "[LOG] " << cc.size() << " connected components (" << smallcomponents << " have less than 100 faces)" << std::endl;
    */

    int cc = tri::Clean<Mesh>::CountConnectedComponents(m);

    Timer t;

    std::unique_ptr<EdgeWeightFunction> wfct(new W3D(m));
    //std::unique_ptr<EdgeWeightFunction> wfct(new W_Geometry3D(m));
    //std::unique_ptr<EdgeWeightFunction> wfct(new WFN(m));
    //std::unique_ptr<EdgeWeightFunction> wfct(new WUV(m));
    std::cout << "[LOG] Weight function " << wfct->Name() << std::endl;
    GraphManager gm{graph, std::move(wfct)};

    double areaThreshold = 0.02;

    std::cout << "Min allowed area percentage while merging = " << areaThreshold << std::endl;

    RecomputeSegmentation(gm, args.regionCount + (cc - 1), areaThreshold);

    int c = ParameterizeGraph(gm, strategy, tolerance);
    if (c > 0) std::cout << "WARNING: " << c << " regions were not parameterized correctly" << std::endl;

    PackingOptions opts = { RasterizationBasedPacker::Parameters::CostFuncEnum::MinWastedSpace, true, true, true, false };
    Pack(gm.Graph(), opts);

    std::cout << "Rendering texture..." << std::endl;
    TextureObjectHandle newTexture = RenderTexture(m, textureObject, args.filter, InterpolationMode::Linear, nullptr);

    graph->textureObject = newTexture;

    PrintParameterizationInfo(graph);
    LogAggregateStats("stats_output", graph, newTexture);
    std::vector<RasterizedParameterizationStats> after = GetRasterizationStats(m, newTexture);
    std::cout << "[LOG] Raster stats after processing" << std::endl;
    LogParameterizationStats(graph, after);


    //ScaleTextureCoordinatesToImage(m, newTexture);
    //LogDistortionStats(graph);

    //ScaleTextureCoordinatesToParameterArea(m, newTexture);

    //GenerateDistortionTextures(m, textureObject);

    std::cout << "Processing took " << t.TimeElapsed() << " seconds" << std::endl;

    std::string savename = "out_" + m.name;
    if (SaveMesh(m, savename.c_str(), newTexture, true) == false) {
        std::cout << "Model not saved correctly" << std::endl;
    }

    GLTerminate();

    return 0;
}

int MainGui(Mesh& m, GraphHandle graph, TextureObjectHandle textureObject, Args args)
{
    GLInit();
    MeshViewer viewer(graph, args);
    viewer.Run();
    GLTerminate();
    return 0;
}

int main(int argc, char *argv[])
{
    // Make sure the executable directory is added to Qt's library path
    {
        QFileInfo executableInfo(argv[0]);
        QCoreApplication::addLibraryPath(executableInfo.dir().absolutePath());
    }

    Args args = parse_args(argc, argv);

    Mesh m;
    TextureObjectHandle textureObject;
    int loadMask;

    if (LoadMesh(m, args.filename.c_str(), textureObject, loadMask) == false) {
        std::cout << "Failed to open mesh" << std::endl;
        std::exit(-1);
    }

    ensure_condition(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);


    GLInit();

    // Log original info

    {
        auto dummyGraph = ComputeParameterizationGraph(m, textureObject);
        PrintParameterizationInfo(dummyGraph);
        LogAggregateStats("stats_input", dummyGraph, textureObject);
        std::vector<RasterizedParameterizationStats> before = GetRasterizationStats(m, textureObject);
        std::cout << "[LOG] Raster stats before processing" << std::endl;
        LogParameterizationStats(dummyGraph, before);
    }

    // Preliminary cleaning

    std::cout << "Cleaning mesh..." << std::endl;

    int dupVert = tri::Clean<Mesh>::RemoveDuplicateVertex(m);
    int zeroArea = tri::Clean<Mesh>::RemoveZeroAreaFace(m);

    tri::Allocator<Mesh>::CompactEveryVector(m);

    tri::UpdateTopology<Mesh>::FaceFace(m);

    if (dupVert > 0)
        std::cout << "Removed " << dupVert << " duplicate vertices" << std::endl;

    if (zeroArea > 0)
        std::cout << "Removed " << zeroArea << " zero area faces" << std::endl;

    int numVertexSplit = tri::Clean<Mesh>::SplitNonManifoldVertex(m, 0);
    if (numVertexSplit > 0)
        std::cout << "Mesh was not vertex manifold, split " << numVertexSplit << " vertices" << std::endl;

    int numRemovedFaces = tri::Clean<Mesh>::RemoveNonManifoldFace(m);
    if (numRemovedFaces > 0)
        std::cout << "Mesh was not edge manifold, removed " << numRemovedFaces << " faces" << std::endl;

    tri::Allocator<Mesh>::CompactEveryVector(m);

    {
        auto dummyGraph = ComputeParameterizationGraph(m, textureObject);

        // Apply topological filter if requested

        if (args.fixContextCapture) {
            std::cout << "WARNING: option --fixcontextcapture detected, attempting to remove topological noise" << std::endl;
            CleanSmallComponents(m, dummyGraph, textureObject, 1e-4);
        }

        // WARNING, THIS NEXT FUNCTION INVALIDATES THE GRAPH

        int outliers = RemoveOutliers(dummyGraph);
        if (outliers)
            std::cout << "Removed " << outliers << " outlier charts" << std::endl;

        dummyGraph = nullptr;
    }

    // Prepare the mesh for processing

    ScaleTextureCoordinatesToImage(m, textureObject);
    ComputeParameterizationScaleInfo(m);
    MarkSeamsAsFaux(m);

    auto graph = ComputeParameterizationGraph(m, textureObject);
    ParameterizeZeroAreaRegions(m, graph);
    StoreWedgeTexCoordAsAttribute(m);

    if (args.gui)
        return MainGui(m, graph, textureObject, args);
    else
        return MainCmd(m, graph, textureObject, args);
}
