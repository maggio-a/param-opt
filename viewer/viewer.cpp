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

    int cc = tri::Clean<Mesh>::CountConnectedComponents(m);

    Timer t;

    std::unique_ptr<EdgeWeightFunction> wfct(new W3D(m));
    LOG_DEBUG << "Weight function " << wfct->Name();
    GraphManager gm{graph, std::move(wfct)};

    double areaThreshold = 0.02;

    LOG_DEBUG << "Min allowed area percentage while merging = " << areaThreshold;

    RecomputeSegmentation(gm, args.regionCount + (cc - 1), areaThreshold);

    int c = ParameterizeGraph(gm, strategy, tolerance);
    if (c > 0)
        LOG_WARN << c << " Regions were not parameterized correctly";

    if (gm.Graph()->Count() < 500) {
        LOG_INFO << "Packing " << gm.Graph()->Count() << " regions in low resolution with permutations";
        PackingOptions opts = { RasterizationBasedPacker::Parameters::CostFuncEnum::MinWastedSpace, true, true, true, false };
        Pack(gm.Graph(), opts);
    } else {
        LOG_INFO << "Packing " << gm.Graph()->Count() << " regions in high resolution";
        PackingOptions opts = { RasterizationBasedPacker::Parameters::CostFuncEnum::MinWastedSpace, false, false, true, false };
        Pack(gm.Graph(), opts);
    }

    LOG_INFO << "Rendering texture...";
    TextureObjectHandle newTexture = RenderTexture(m, textureObject, args.filter, InterpolationMode::Linear, nullptr);

    graph->textureObject = newTexture;

    PrintParameterizationInfo(graph);
    LogAggregateStats("stats_output", graph, newTexture);
    std::vector<RasterizedParameterizationStats> after = GetRasterizationStats(m, newTexture);
    LOG_INFO << "Raster stats after processing";
    LogParameterizationStats(graph, after);

    LOG_VERBOSE << "Processing took " << t.TimeElapsed() << " seconds";

    std::string savename = "out_" + m.name;
    if (SaveMesh(m, savename.c_str(), newTexture, true) == false) {
        LOG_ERR << "Model not saved correctly";
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

    LOG_INIT(args.logLevel);

    Mesh m;
    TextureObjectHandle textureObject;
    int loadMask;

    if (LoadMesh(m, args.filename.c_str(), textureObject, loadMask) == false) {
        LOG_ERR << "Failed to open mesh";
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
        LOG_INFO << "Raster stats before processing";
        LogParameterizationStats(dummyGraph, before);
    }

    // Preliminary cleaning

    LOG_INFO << "Cleaning mesh...";

    int dupVert = tri::Clean<Mesh>::RemoveDuplicateVertex(m);
    int zeroArea = tri::Clean<Mesh>::RemoveZeroAreaFace(m);

    tri::Allocator<Mesh>::CompactEveryVector(m);

    tri::UpdateTopology<Mesh>::FaceFace(m);

    if (dupVert > 0)
        LOG_VERBOSE << "Removed " << dupVert << " duplicate vertices";

    if (zeroArea > 0)
        LOG_VERBOSE << "Removed " << zeroArea << " zero area faces";

    int numVertexSplit = 0;
    int nv;
    while ((nv = tri::Clean<Mesh>::SplitNonManifoldVertex(m, 0)) > 0)
        numVertexSplit += nv;
    if (numVertexSplit > 0)
        LOG_VERBOSE << "Mesh was not vertex manifold, split " << numVertexSplit << " vertices";

    int numRemovedFaces = tri::Clean<Mesh>::RemoveNonManifoldFace(m);
    if (numRemovedFaces > 0)
        LOG_VERBOSE << "Mesh was not edge manifold, removed " << numRemovedFaces << " faces";

    tri::Allocator<Mesh>::CompactEveryVector(m);

    {
        auto dummyGraph = ComputeParameterizationGraph(m, textureObject);

        // Apply topological filter if requested

        if (args.fixContextCapture) {
            LOG_WARN << "Option --fixcontextcapture detected, attempting to remove topological noise";
            CleanSmallComponents(m, dummyGraph, textureObject, 1e-4);
        }

        // WARNING, THIS NEXT FUNCTION INVALIDATES THE GRAPH

        int outliers = RemoveOutliers(dummyGraph);
        if (outliers)
            LOG_VERBOSE << "Removed " << outliers << " outlier charts";

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
