#include "mesh.h"
#include "mesh_graph.h"
#include "uv.h"
#include "texture_optimization.h"
#include "texture_rendering.h"
#include "parameterization_checker.h"
#include "timer.h"
#include "gl_utils.h"
#include "texture.h"
#include "mesh_viewer.h"
#include "logging.h"

#include <wrap/io_trimesh/io_mask.h>

#include <string>
#include <vector>
#include <iostream>

#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>


#define CLEAN_TROUBLESOME_REGIONS_CONTEXTCAPTURE


using namespace vcg;

struct Args {
    std::string filename;
    bool gui;
    bool filter;

    Args() : filename{}, gui{false}, filter{true} {}
};

Args parse_args(int argc, char *argv[])
{
    Args args;
    for (int i = 1; i < argc; ++i) {
        std::string s(argv[i]);
        if (s.substr(0, 2) == std::string("--")) {
            if (s == std::string("--gui"))
                args.gui = true;
            else if (s == std::string("--nofilter"))
                args.filter = false;
            else
                assert(0 && "Invalid flag option");
        } else {
            assert(args.filename.empty() && "Invalid argument");
            args.filename = s;
        }
    }
    return args;
}

int MainCmd(Mesh& m, GraphHandle graph, TextureObjectHandle textureObject,
             Args args)
{
    int minRegionSize;
    if (m.FN() < 100000)
        minRegionSize = 1000;
    else if (m.FN() < 300000)
        minRegionSize = 5000;
    else
        minRegionSize = 10000;

    if (minRegionSize > m.FN()) {
        std::cout << "WARNING: minFaceCount > m.FN()" << std::endl;
    }

    ParameterizationStrategy strategy = MakeStrategy(
            DirectParameterizer::FixedBorderBijective,
            EnergyType::SymmetricDirichlet,
            ParameterizationGeometry::Texture,
            //DescentType::ScalableLocallyInjectiveMappings,
            DescentType::CompositeMajorization,
            500,            // Number of iterations
            true,           // Fill holes ?
            true,           // Use cuts ?
            false,          // Use warm start ?
            true            // Use scaffolding ?
    );
    double tolerance = 0.0002;

    LogStrategy(strategy, tolerance);

    std::vector<std::pair<int, Mesh::FacePointer>> cc;
    tri::Clean<Mesh>::ConnectedComponents(m, cc);
    int smallcomponents = 0;
    for (auto p : cc) if (p.first < 100) smallcomponents++;

    std::cout << "[LOG] " << cc.size() << " connected components (" << smallcomponents << " have less than 100 faces)" << std::endl;

    Timer t;

    std::unique_ptr<EdgeWeightFunction> wfct(new W3D(m));
    //std::unique_ptr<EdgeWeightFunction> wfct(new WFN(m));
    //std::unique_ptr<EdgeWeightFunction> wfct(new WUV(m));
    std::cout << "[LOG] Weight function " << wfct->Name() << std::endl;
    GraphManager gm{graph, std::move(wfct)};

    int regionCount = 20;
    RecomputeSegmentation(gm, regionCount + smallcomponents, minRegionSize);

    int c = ParameterizeGraph(gm, strategy, tolerance);
    if (c > 0) std::cout << "WARNING: " << c << " regions were not parameterized correctly" << std::endl;

    PackingOptions opts = { RasterizationBasedPacker::Parameters::CostFuncEnum::MinWastedSpace, true, true, true, false };
    Pack(gm.Graph(), opts);

    std::cout << "Rendering texture..." << std::endl;
    TextureObjectHandle newTexture = RenderTexture(m, textureObject, args.filter, InterpolationMode::Linear, nullptr);

    std::vector<RasterizedParameterizationStats> after = GetRasterizationStats(m, newTexture);
    std::cout << "[LOG] Raster stats after processing" << std::endl;
    LogParameterizationStats(graph, after);


    ScaleTextureCoordinatesToImage(m, newTexture);
    LogDistortionStats(graph);

    ScaleTextureCoordinatesToParameterArea(m, newTexture);

    //GenerateDistortionTextures(m, textureObject);

    std::cout << "Processing took " << t.TimeElapsed() << " seconds" << std::endl;

    std::string savename = "out_" + m.name;
    if (SaveMesh(m, savename.c_str(), newTexture, true) == false) {
        std::cout << "Model not saved correctly" << std::endl;
    }

    PrintParameterizationInfo(graph);


    GLTerminate();

    return 0;
}

int MainGui(Mesh& m, GraphHandle graph, TextureObjectHandle textureObject,
             Args args)
{
    GLInit();
    MeshViewer viewer(graph, args.filename);
    viewer.Run();
    GLTerminate();
    return 0;
}

int main(int argc, char *argv[])
{
    std::cout << "sizeof Mesh = " << sizeof(Mesh) << std::endl;
    std::cout << "sizeof MeshFace = " << sizeof(MeshFace) << std::endl;
    std::cout << "sizeof MeshVertex = " << sizeof(MeshVertex) << std::endl;

    Args args = parse_args(argc, argv);

    Mesh m;
    TextureObjectHandle textureObject;
    int loadMask;

    if (LoadMesh(m, args.filename.c_str(), textureObject, loadMask) == false) {
        std::cout << "Failed to open mesh" << std::endl;
        std::exit(-1);
    }

    assert(textureObject->ArraySize() == 1 && "Currently only single texture is supported");
    assert(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);


    // Print original info
    auto dummyGraph = ComputeParameterizationGraph(m, textureObject);
    PrintParameterizationInfo(dummyGraph);

    GLInit();
    std::vector<RasterizedParameterizationStats> before = GetRasterizationStats(m, textureObject);
    std::cout << "[LOG] Raster stats before processing" << std::endl;
    LogParameterizationStats(dummyGraph, before);


    tri::UpdateTopology<Mesh>::FaceFace(m);
    int numVertexSplit = tri::Clean<Mesh>::SplitNonManifoldVertex(m, 0);
    if (numVertexSplit > 0)
        std::cout << "[LOG] Mesh was not vertex manifold, split " << numVertexSplit << " vertices" << std::endl;
    int numRemovedFaces = tri::Clean<Mesh>::RemoveNonManifoldFace(m);
    if (numRemovedFaces > 0)
        std::cout << "[LOG] Mesh was not edge manifold, removed " << numRemovedFaces << " faces" << std::endl;

    int smallCCSizeThreshold = int(2.5e-4 * m.FN());
    auto p = tri::Clean<Mesh>::RemoveSmallConnectedComponentsSize(m, 2.5e-4 * m.FN());
    std::cout << "[LOG] RemoveSmallConnectedComponents (total cc, removed cc, threshold) = "
              << p.first << " , "  << p.second << " , " << smallCCSizeThreshold << std::endl;

    tri::Allocator<Mesh>::CompactEveryVector(m);


#ifdef CLEAN_TROUBLESOME_REGIONS_CONTEXTCAPTURE
    CleanSmallComponents(m, dummyGraph, textureObject, 1e-4);
#endif

    ScaleTextureCoordinatesToImage(m, textureObject);

    dummyGraph = nullptr;

    ComputeParameterizationScaleInfo(m);
    MarkSeamsAsFaux(m);

    auto graph = ComputeParameterizationGraph(m, textureObject);
    PreprocessMesh(m, graph);
    StoreWedgeTexCoordAsAttribute(m);


    // Print original info
    PrintParameterizationInfo(graph);

    if (args.gui)
        return MainGui(m, graph, textureObject, args);
    else
        return MainCmd(m, graph, textureObject, args);
}
