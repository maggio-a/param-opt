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

#include <wrap/io_trimesh/io_mask.h>

#include <string>
#include <vector>
#include <iostream>

#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>


//#define CLEAN_TROUBLESOME_REGIONS_CONTEXTCAPTURE


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


void LogStrategy(ParameterizationStrategy strategy, double tol)
{
    std::string geometry = (strategy.geometry == ParameterizationGeometry::Model) ? std::string("Model") : std::string("Texture");
    std::cout << "[LOG] Parameterization strategy: ";
    std::cout << "Geometry=" << geometry << " , "
              << "iterations=" << strategy.optimizerIterations << " , "
              << "padded inner boundaries=" << strategy.padBoundaries << " , "
              << "cuts=" << strategy.applyCut << " , "
              << "scaffold=" << strategy.scaffold << " , "
              << "tolerance=" << tol
              << std::endl;
}

void PercentilePlot(vcg::Distribution<double>& d)
{
    constexpr int PERC_COUNT = 20;
    for (int i = 0; i < PERC_COUNT; ++i) {
        std::cout << "[LOG], ";
        std::cout << (1.0 / 20) * i << " , ";
        std::cout << d.Percentile((1.0 / 20) * i) << std::endl;
    }
    std::cout << "[LOG], ";
    std::cout << 1.0 << " , ";
    std::cout << d.Percentile(1.0) << std::endl;
}

void LogDistortionStats(std::shared_ptr<MeshGraph> graph, double uvRatio)
{
    assert(0 && "multitexture not fixed");
    //  Correct ratio --------------------------------------------------
    if (uvRatio != 1.0) {
        auto wtcsattr = GetWedgeTexCoordStorageAttribute(graph->mesh);
        for (auto& f : graph->mesh.face) {
            for (int i = 0; i < 3; ++i) {
                wtcsattr[f].tc[i].U() /= (uvRatio > 1.0) ? uvRatio : 1.0;
                wtcsattr[f].tc[i].V() *= (uvRatio < 1.0) ? uvRatio : 1.0;
            }
        }
    }

    vcg::Distribution<double> d;

    std::cout << "[LOG] Distortion" << std::endl;
    std::cout << "[LOG] Type Source Min PCT1 PCT5 Max PCT99 PCT95 Avg Variance" << std::endl;

    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Area, ParameterizationGeometry::Texture);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Area Texture "
              << d.Min() << " " << d.Percentile(0.01) << " " << d.Percentile(0.05) << " "
              << d.Max() << " " << d.Percentile(0.99) << " " << d.Percentile(0.95) << " "
              << d.Avg() << " " << d.Variance() << std::endl;

    std::cout << "[LOG]  Percentile plot: " << std::endl;
    PercentilePlot(d);

    d.Clear();

    graph->MapDistortion(DistortionMetric::Type::Angle, ParameterizationGeometry::Texture);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Angle Texture "
              << d.Min() << " " << d.Percentile(0.01) << " " << d.Percentile(0.05) << " "
              << d.Max() << " " << d.Percentile(0.99) << " " << d.Percentile(0.95) << " "
              << d.Avg() << " " << d.Variance() << std::endl;
    std::cout << "[LOG]  Percentile plot: " << std::endl;
    PercentilePlot(d);

    //  Correct ratio --------------------------------------------------
    if (uvRatio != 1.0) {
        auto wtcsattr = GetWedgeTexCoordStorageAttribute(graph->mesh);
        for (auto& f : graph->mesh.face) {
            for (int i = 0; i < 3; ++i) {
                wtcsattr[f].tc[i].U() *= (uvRatio > 1.0) ? uvRatio : 1.0;
                wtcsattr[f].tc[i].V() /= (uvRatio < 1.0) ? uvRatio : 1.0;
            }
        }
    }

    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Area, ParameterizationGeometry::Model);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Area Model "
              << d.Min() << " " << d.Percentile(0.01) << " " << d.Percentile(0.05) << " "
              << d.Max() << " " << d.Percentile(0.99) << " " << d.Percentile(0.95) << " "
              << d.Avg() << " " << d.Variance() << std::endl;
    std::cout << "[LOG]  Percentile plot: " << std::endl;
    PercentilePlot(d);

    d.Clear();

    if (uvRatio != 1.0) {
        for (auto& f : graph->mesh.face) {
            for (int i = 0; i < 3; ++i) {
                f.cWT(i).U() *= (uvRatio > 1.0) ? uvRatio : 1.0;
                f.cWT(i).V() /= (uvRatio < 1.0) ? uvRatio : 1.0;
            }
        }
    }

    graph->MapDistortion(DistortionMetric::Type::Angle, ParameterizationGeometry::Model);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Angle Model "
              << d.Min() << " " << d.Percentile(0.01) << " " << d.Percentile(0.05) << " "
              << d.Max() << " " << d.Percentile(0.99) << " " << d.Percentile(0.95) << " "
              << d.Avg() << " " << d.Variance() << std::endl;
    std::cout << "[LOG]  Percentile plot: " << std::endl;
    PercentilePlot(d);

    if (uvRatio != 1.0) {
        for (auto& f : graph->mesh.face) {
            for (int i = 0; i < 3; ++i) {
                f.cWT(i).U() /= (uvRatio > 1.0) ? uvRatio : 1.0;
                f.cWT(i).V() *= (uvRatio < 1.0) ? uvRatio : 1.0;
            }
        }
    }
}

void LogParameterizationStats(std::shared_ptr<MeshGraph> graph, RasterizedParameterizationStats stats, const std::string& header)
{
    int boundaries;
    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(graph->mesh);
    boundaries = tri::Clean<Mesh>::CountHoles(graph->mesh);
    tri::UpdateTopology<Mesh>::FaceFace(graph->mesh);

    std::cout << header << std::endl;
    std::cout << "[LOG] Texture Size , Charts , Boundaries , Occupancy , BorderUV(px), AreaUV(px), BilinearAreaUV(px), Overwritten fragments, Lost Fragments" << std::endl;
    std::cout << "[LOG] "
              << stats.rw << "x" << stats.rh << " , "
              << graph->charts.size() << " , "
              << boundaries << " , "
              << (stats.totalFragments - stats.lostFragments) / double(stats.rw*stats.rh) << " , "
              << stats.boundaryFragments << " , "
              << stats.totalFragments  << " , "
              << stats.totalFragments_bilinear  << " , "
              << stats.overwrittenFragments  << " , "
              << stats.lostFragments << "" << std::endl;

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
            DescentType::ScalableLocallyInjectiveMappings,
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

    GLInit();

    RasterizedParameterizationStats before = GetRasterizationStats(m, textureObject->imgVec[0]->width(), textureObject->imgVec[0]->height());
    LogParameterizationStats(graph, before, std::string("[LOG] Raster stats before parameterizing"));

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

    RasterizedParameterizationStats after = GetRasterizationStats(m, textureObject->imgVec[0]->width(), textureObject->imgVec[0]->height());
    LogParameterizationStats(graph, after, std::string("[LOG] Raster stats after parameterizing"));

    std::cout << "Rendering texture..." << std::endl;
    TextureObjectHandle newTexture = RenderTexture(m, textureObject, args.filter, InterpolationMode::Linear, nullptr);

    //GenerateDistortionTextures(m, textureObject);

    std::cout << "Processing took " << t.TimeElapsed() << " seconds" << std::endl;

    std::string savename = "out_" + m.name;
    if (SaveMesh(m, savename.c_str(), newTexture, true) == false) {
        std::cout << "Model not saved correctly" << std::endl;
    }

    PrintParameterizationInfo(graph);

    //LogDistortionStats(graph, newTexture->TextureWidth(0) / (double) newTexture->TextureHeight(0));

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

    // FIXME this seems to slow down the parameterization process... conditioning issues?
    ScaleTextureCoordinatesToImage(m, textureObject);

    // Print original info
    auto dummyGraph = ComputeParameterizationGraph(m, textureObject);
    PrintParameterizationInfo(dummyGraph);

#ifdef CLEAN_TROUBLESOME_REGIONS_CONTEXTCAPTURE
    CleanSmallComponents(m, dummyGraph, textureObject, 1e-4);
#endif

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
