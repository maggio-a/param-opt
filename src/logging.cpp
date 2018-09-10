#include "logging.h"
#include "parameterization.h"
#include "texture_rendering.h"

void LogStrategy(const ParameterizationStrategy& strategy, double tol)
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
}

void LogParameterizationStats(std::shared_ptr<MeshGraph> graph, const std::vector<RasterizedParameterizationStats>& statsVec)
{
    int boundaries;
    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(graph->mesh);
    boundaries = tri::Clean<Mesh>::CountHoles(graph->mesh);
    tri::UpdateTopology<Mesh>::FaceFace(graph->mesh);

    long rasterArea = 0;
    long usedFragments = 0;
    long overwrittenFragments = 0;
    long lostFragments = 0;
    long boundaryFragments = 0;

    for (const auto& stats : statsVec) {
        rasterArea += (stats.rw * stats.rh);
        usedFragments += (stats.totalFragments - stats.lostFragments);
        boundaryFragments += stats.boundaryFragments;
        overwrittenFragments += stats.overwrittenFragments;
        lostFragments += stats.lostFragments;
    }

    std::cout << "[LOG] Texture Area , Charts , Boundaries , Occupancy , Border, Lost Fragments" << std::endl;
    std::cout << "[LOG] "
              << rasterArea << " , "
              << graph->Count() << " , "
              << boundaries << " , "
              << usedFragments / double(rasterArea) << " , "
              << graph->BorderUV() << " , "
              << lostFragments << "" << std::endl;

}

