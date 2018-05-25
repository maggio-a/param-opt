#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/io_mask.h>

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
#include "parameterization_checker.h"
#include "timer.h"
#include "gl_utils.h"
#include "texture.h"

#include "mesh_viewer.h"

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
              << "tolerance=" << tol << " , "
              << "padded inner boundaries=" << strategy.padBoundaries
              << std::endl;
}

void LogDistortionStats(std::shared_ptr<MeshGraph> graph)
{
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

    std::cout << "[LOG]  Percentile plot: ";
    for (double p = 0.0; p <= 1.0; p += 0.05) {
        std::cout << "( " << p << " , " << d.Percentile(p) << " ) ; ";
    }
    std::cout << "( " << 1.0 << " , " << d.Percentile(1.0) << " ) ; ";
    std::cout << endl;

    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Angle, ParameterizationGeometry::Texture);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Angle Texture "
              << d.Min() << " " << d.Percentile(0.01) << " " << d.Percentile(0.05) << " "
              << d.Max() << " " << d.Percentile(0.99) << " " << d.Percentile(0.95) << " "
              << d.Avg() << " " << d.Variance() << std::endl;
    std::cout << "[LOG]  Percentile plot: ";
    for (double p = 0.0; p <= 1.0; p += 0.05) {
        std::cout << "( " << p << " , " << d.Percentile(p) << " ) ; ";
    }
    std::cout << "( " << 1.0 << " , " << d.Percentile(1.0) << " ) ; ";
    std::cout << endl;


    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Area, ParameterizationGeometry::Model);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Area Model "
              << d.Min() << " " << d.Percentile(0.01) << " " << d.Percentile(0.05) << " "
              << d.Max() << " " << d.Percentile(0.99) << " " << d.Percentile(0.95) << " "
              << d.Avg() << " " << d.Variance() << std::endl;
    std::cout << "[LOG]  Percentile plot: ";
    for (double p = 0.0; p <= 1.0; p += 0.05) {
        std::cout << "( " << p << " , " << d.Percentile(p) << " ) ; ";
    }
    std::cout << "( " << 1.0 << " , " << d.Percentile(1.0) << " ) ; ";
    std::cout << endl;

    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Angle, ParameterizationGeometry::Model);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Angle Model "
              << d.Min() << " " << d.Percentile(0.01) << " " << d.Percentile(0.05) << " "
              << d.Max() << " " << d.Percentile(0.99) << " " << d.Percentile(0.95) << " "
              << d.Avg() << " " << d.Variance() << std::endl;
    std::cout << "[LOG]  Percentile plot: ";
    for (double p = 0.0; p <= 1.0; p += 0.05) {
        std::cout << "( " << p << " , " << d.Percentile(p) << " ) ; ";
    }
    std::cout << "( " << 1.0 << " , " << d.Percentile(1.0) << " ) ; ";
    std::cout << endl;

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

int main_cmd(Mesh& m, GraphHandle graph, TextureObjectHandle textureObject,
             Args args)
{
    int minRegionSize;
    if (m.FN() < 100000)
        minRegionSize = 1000;
    else if (m.FN() < 300000)
        minRegionSize = 5000;
    else
        minRegionSize = 10000;

    //assert(m.FN() < 2000000);
    //assert(m.FN() < 1000000);

    if (minRegionSize > m.FN()) {
        std::cout << "WARNING: minFaceCount > m.FN()" << std::endl;
    }

    ParameterizationStrategy strategy;
    strategy.directParameterizer = FixedBorderBijective;
    strategy.energy = EnergyType::SymmetricDirichlet;
    strategy.geometry = Texture;
    //strategy.geometry = Model;
    strategy.descent = ScalableLocallyInjectiveMappings;
    strategy.optimizerIterations = 500;
    strategy.padBoundaries = true;
    strategy.applyCut = true;
    strategy.warmStart = false;
    double tolerance = 0.0005;

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
    ReduceTextureFragmentation_NoPacking_TargetRegionCount(gm, regionCount + smallcomponents, minRegionSize);
    //ReduceTextureFragmentation_NoPacking(gm, minRegionSize);

    int c = ParameterizeGraph(gm, 1.0, strategy, true, tolerance, true);
    if (c > 0) std::cout << "WARNING: " << c << " regions were not parameterized correctly" << std::endl;

    LogDistortionStats(graph);

    RasterizedParameterizationStats after = GetRasterizationStats(m, textureObject->imgVec[0]->width(), textureObject->imgVec[0]->height());
    LogParameterizationStats(graph, after, std::string("[LOG] Raster stats after parameterizing"));

    std::cout << "Rendering texture..." << std::endl;
    TextureObjectHandle newTexture = RenderTexture(m, textureObject, args.filter, true, nullptr);

    std::cout << "Processing took " << t.TimeElapsed() << " seconds" << std::endl;

    std::string savename = "out_" + m.name;
    if (SaveMesh(m, savename.c_str(), newTexture, true) == false) {
        std::cout << "Model not saved correctly" << std::endl;
    }

    PrintParameterizationInfo(graph);

    GLTerminate();

    return 0;
}

int main_gui(Mesh& m, GraphHandle graph, TextureObjectHandle textureObject,
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

    //double packingCoverage;
    //CompactTextureData(m, textureObject, &packingCoverage);

    tri::UpdateTopology<Mesh>::FaceFace(m);
    if (tri::Clean<Mesh>::CountNonManifoldEdgeFF(m, false)) {
        std::cout << "Mesh is not edge manifold" << std::endl;
        std::exit(-1);
    }

    ComputeParameterizationScaleInfo(m);
    MarkSeamsAsFaux(m);
    StoreWedgeTexCoordAsAttribute(m);
    PreprocessMesh(m);

    auto graph = ComputeParameterizationGraph(m, textureObject);

    // Print original info
    PrintParameterizationInfo(graph);

    if (args.gui)
        return main_gui(m, graph, textureObject, args);
    else
        return main_cmd(m, graph, textureObject, args);
}
