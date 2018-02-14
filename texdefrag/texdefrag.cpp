#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/export.h>

#include <string>
#include <vector>
#include <iostream>
#include <memory>

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

using namespace vcg;

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


    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Angle, ParameterizationGeometry::Texture);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Angle Texture "
              << d.Min() << " " << d.Percentile(0.01) << " " << d.Percentile(0.05) << " "
              << d.Max() << " " << d.Percentile(0.99) << " " << d.Percentile(0.95) << " "
              << d.Avg() << " " << d.Variance() << std::endl;

    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Area, ParameterizationGeometry::Model);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Area Model "
              << d.Min() << " " << d.Percentile(0.01) << " " << d.Percentile(0.05) << " "
              << d.Max() << " " << d.Percentile(0.99) << " " << d.Percentile(0.95) << " "
              << d.Avg() << " " << d.Variance() << std::endl;

    d.Clear();
    graph->MapDistortion(DistortionMetric::Type::Angle, ParameterizationGeometry::Model);
    for (auto& f : graph->mesh.face) {
        d.Add(f.Q());
    }
    std::cout << "[LOG] Angle Model "
              << d.Min() << " " << d.Percentile(0.01) << " " << d.Percentile(0.05) << " "
              << d.Max() << " " << d.Percentile(0.99) << " " << d.Percentile(0.95) << " "
              << d.Avg() << " " << d.Variance() << std::endl;
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
    //strategy.geometry = Model;
    strategy.descent = ScalableLocallyInjectiveMappings;
    strategy.optimizerIterations = 300;
    strategy.padBoundaries = true;

    double tolerance = 0.0005;

    double initialArea = 0.0;
    std::unordered_set<std::size_t> zeroIndices;
    for (auto& f : m.face) {
        double area = DistortionMetric::AreaUV(f);
        if (area == 0) {
            zeroIndices.insert(tri::Index(m, f));
        }
        else initialArea += area;
    }

    PreprocessMesh(m);
    StoreWedgeTexCoordAsAttribute(m);

    LogStrategy(strategy, tolerance);

    std::vector<std::pair<int, Mesh::FacePointer>> cc;
    tri::Clean<Mesh>::ConnectedComponents(m, cc);
    int smallcomponents = 0;
    for (auto p : cc) if (p.first < 100) smallcomponents++;

    std::cout << "[LOG] " << cc.size() << " connected components (" << smallcomponents << " have less than 100 faces)" << std::endl;

    float uvMeshBorder;
    auto graph = ComputeParameterizationGraph(m, textureObject, &uvMeshBorder);

    double areaUvBefore = graph->AreaUV();

    GLInit();

    // Print original info
    PrintParameterizationInfo(graph);

    RasterizedParameterizationStats before = GetRasterizationStats(m, textureObject->imgVec[0]->width(), textureObject->imgVec[0]->height());
    LogParameterizationStats(graph, before, std::string("[LOG] Raster stats before parameterizing"));

    Timer t;

    std::unique_ptr<EdgeWeightFunction> wfct(new FaceSizeWeightedShared3DBorder(m));
    std::cout << "[LOG] Weight function" << wfct->Name() << std::endl;
    GraphManager gm{graph, std::move(wfct)};

    int regionCount = 20;
    if (regionCount > 0) {
        ReduceTextureFragmentation_NoPacking_TargetRegionCount(gm, regionCount + smallcomponents, minRegionSize);
    }
    else {
        ReduceTextureFragmentation_NoPacking(gm, minRegionSize);
    }

    int c = ParameterizeGraph(gm, strategy, true, tolerance, true);
    if (c > 0) std::cout << "WARNING: " << c << " regions were not parameterized correctly" << std::endl;

    double finalArea = 0.0;
    for (auto& f : m.face) {
        double area = DistortionMetric::AreaUV(f);
        if (zeroIndices.count(tri::Index(m, f)) == 0) {
            finalArea += area;
        }
    }


    std::cout << "[LOG] Scale factor of the packed chart = " << initialArea / finalArea << std::endl;
    //std::cout << "[LOG] Scale factor of the packed chart = " << graph->AreaUV() / areaUvBefore << std::endl;

    bool normalizeArea = true;
    //scale parameterization
    double scale = std::sqrt(initialArea / finalArea);
    if (normalizeArea)  {
        std::cout << "[LOG] *** AREA NORMALIZED *** " << std::endl;
        for (auto& f : m.face) {
            f.WT(0).P() *= scale;
            f.WT(1).P() *= scale;
            f.WT(2).P() *= scale;
        }
    }

    LogDistortionStats(graph);

    RasterizedParameterizationStats after = GetRasterizationStats(m, textureObject->imgVec[0]->width(), textureObject->imgVec[0]->height());
    LogParameterizationStats(graph, after, std::string("[LOG] Raster stats after parameterizing"));

    std::cout << "Rendering texture..." << std::endl;
    TextureObjectHandle newTexture = RenderTexture(m, textureObject, filter, nullptr);

    std::cout << "Processing took " << t.TimeElapsed() << " seconds" << std::endl;

    graph->MapDistortion(DistortionMetric::Type::Angle, ParameterizationGeometry::Texture);

    std::string outName = "out_" + modelName;
    if (SaveMesh(m, outName.c_str(), newTexture, true) == false) {
        std::cout << "Model not saved correctly" << std::endl;
    }

    auto graph2 = ComputeParameterizationGraph(m, textureObject, &uvMeshBorder);
    // Print optimized info
    PrintParameterizationInfo(graph2);

    GLTerminate();

    return 0;
}


