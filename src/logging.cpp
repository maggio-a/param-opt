#include "logging.h"
#include "parameterization.h"
#include "texture_rendering.h"
#include "utils.h"

#include <Eigen/Core>

#include <vcg/math/histogram.h>

#include <string>
#include <fstream>

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

void ExtractSingularValues(vcg::Point3d p10, vcg::Point3d p20, vcg::Point2d u10, vcg::Point2d u20, double *s_min, double *s_max)
{
    Point2d x10;
    Point2d x20;
    LocalIsometry(p10, p20, x10, x20);

    // compute the singular values of the transformation matrix, s2 > s1
    // ref for the formula: smith&schaefer 2015 bijective,
    Eigen::Matrix2d phi = ComputeTransformationMatrix(x10, x20, u10, u20);
    double bcplus  = std::pow(phi(0, 1) + phi(1, 0), 2.0);
    double bcminus = std::pow(phi(0, 1) - phi(1, 0), 2.0);
    double adplus  = std::pow(phi(0, 0) + phi(1, 1), 2.0);
    double adminus = std::pow(phi(0, 0) - phi(1, 1), 2.0);
    *s_min = 0.5 * std::abs(std::sqrt(bcplus + adminus) - std::sqrt(bcminus + adplus));
    *s_max = 0.5 * (std::sqrt(bcplus + adminus) + std::sqrt(bcminus + adplus));
}

static std::string JSONField(const char *fieldName, int fieldValue)
{
    std::stringstream ss;
    ss << "\"" << fieldName << "\": " << fieldValue;
    return ss.str();
}

static std::string JSONField(const char *fieldName, double fieldValue)
{
    std::stringstream ss;
    ss << "\"" << fieldName << "\": " << fieldValue;
    return ss.str();
}

static std::string ToJSON(const RasterizedParameterizationStats& stats)
{
    std::stringstream ss;
    ss << "{" << std::endl;
    ss << JSONField("rw"                         , stats.rw) << "," << std::endl;
    ss << JSONField("rh"                         , stats.rh) << "," << std::endl;
    ss << JSONField("totalFragments"             , stats.totalFragments) << "," << std::endl;
    ss << JSONField("totalFragments_bilinear"    , stats.totalFragments_bilinear) << "," << std::endl;
    ss << JSONField("overwrittenFragments"       , stats.overwrittenFragments) << "," << std::endl;
    ss << JSONField("lostFragments"              , stats.lostFragments) << "," << std::endl;
    ss << JSONField("fragmentClashes"            , stats.fragmentClashes) << "," << std::endl;
    ss << JSONField("boundaryFragments"          , stats.boundaryFragments) << std::endl;
    ss << "}";

    return ss.str();
}

static std::string JSONField(const char *fieldName, const std::vector<RasterizedParameterizationStats>& vs)
{
    std::stringstream ss;
    ss << "\"" << fieldName << "\": [" << std::endl;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i > 0)
            ss << "," << std::endl;
        ss << ToJSON(vs[i]);
    }
    ss << std::endl << "]";

    return ss.str();
}

static std::string ToJSON(const GeometryImageStats& stats)
{
    std::stringstream ss;
    ss << "{" << std::endl;
    ss << JSONField("rw"            , stats.rw) << "," << std::endl;
    ss << JSONField("rh"            , stats.rh) << "," << std::endl;
    ss << JSONField("uniqueFaults"  , stats.uniqueFaults) << "," << std::endl;
    ss << JSONField("totalFaults"   , stats.totalFaults) << std::endl;
    ss << "}";

    return ss.str();
}

static std::string JSONField(const char *fieldName, const std::vector<GeometryImageStats>& gis)
{
    std::stringstream ss;
    ss << "\"" << fieldName << "\": [" << std::endl;
    for (std::size_t i = 0; i < gis.size(); ++i) {
        if (i > 0)
            ss << "," << std::endl;
        ss << ToJSON(gis[i]);
    }
    ss << std::endl << "]";
    return ss.str();
}

static std::string JSONField(const char *fieldName, const char *str)
{
    std::stringstream ss;
    ss << "\"" << fieldName << "\": \"" << str << "\"";
    return ss.str();
}

static std::string JSONField(const char *fieldName, bool b)
{
    std::stringstream ss;
    ss << "\"" << fieldName << "\": " << (b ? "true" : "false");
    return ss.str();
}

static std::string JSONField(const char *fieldName, const std::string& str)
{
    std::stringstream ss;
    ss << "\"" << fieldName << "\": \"" << str << "\"";
    return ss.str();
}

static std::string JSONField(const char *fieldName, const std::vector<double>& v)
{
    std::stringstream ss;
    ss << "\"" << fieldName << "\": [" << std::endl;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i > 0)
            ss << ", ";
        ss << v[i];
    }
    ss << std::endl << "]";
    return ss.str();
}

static std::string ToJSON(const vcg::Histogram<double>& hist)
{
    std::stringstream ss;
    ss << "{" << std::endl;
    ss << JSONField("H", hist.H) << "," << std::endl;
    ss << JSONField("R", hist.R) << "," << std::endl;
    ss << JSONField("minv", hist.minv) << "," << std::endl;
    ss << JSONField("maxv", hist.maxv) << "," << std::endl;
    ss << JSONField("minElem", hist.minElem) << "," << std::endl;
    ss << JSONField("maxElem", hist.maxElem) << "," << std::endl;
    ss << JSONField("n", hist.n) << "," << std::endl;
    ss << JSONField("cnt", hist.cnt) << "," << std::endl;
    ss << JSONField("sum", hist.sum) << "," << std::endl;
    ss << JSONField("rms", hist.rms) << std::endl;
    ss << "}";
    return ss.str();
}

static std::string JSONField(const char *fieldName, const vcg::Histogram<double>& hist)
{
    std::stringstream ss;
    ss << "\"" << fieldName << "\": " << ToJSON(hist);
    return ss.str();
}

void LogAggregateStats(const std::string& filename, std::shared_ptr<MeshGraph> graph, TextureObjectHandle textureObject)
{

    Mesh& m = graph->mesh;

    tri::UpdateTopology<Mesh>::FaceFace(m);

    int nonmanif_vert = tri::Clean<Mesh>::CountNonManifoldVertexFF(m);
    int nonmanif_edge = tri::Clean<Mesh>::CountNonManifoldEdgeFF(m);
    int connected_components = tri::Clean<Mesh>::CountConnectedComponents(m);
    int boundary_loops = tri::Clean<Mesh>::CountHoles(m);
    int genus = tri::Clean<Mesh>::MeshGenus(m);

    double mapped_fraction = graph->MappedFraction();


    // number of chart
    int num_charts = graph->Count();


    // number of null chart
    int num_null_charts = 0;
    for (auto entry : graph->charts)
        if (entry.second->AreaUV() == 0)
            num_null_charts++;


    // angle distortion histogram
    ScaleTextureCoordinatesToImage(m, textureObject);

    vcg::Histogram<double> hist_angle;
    hist_angle.Clear();
    hist_angle.SetRange(1.0, 10.0, 1000);
    for (auto& f : m.face) {
        if ((DistortionMetric::Area3D(f) > 0)) {
            double areaUV = DistortionMetric::AreaUV(f);
            if (std::isfinite(areaUV)) {
                double s_min, s_max;
                ExtractSingularValues(f.P(1) - f.P(0), f.P(2) - f.P(0), f.WT(1).P() - f.WT(0).P(), f.WT(2).P() - f.WT(0).P(), &s_min, &s_max);
                double quasi_conformal_distortion = s_max / s_min;
                if (std::isfinite(quasi_conformal_distortion))
                    hist_angle.Add(quasi_conformal_distortion, DistortionMetric::Area3D(f));
            }
        }
    }


    // rasterization stats at various mipmap levels
    ScaleTextureCoordinatesToParameterArea(m, textureObject);

    vcg::Histogram<double> hist_area;
    hist_area.Clear();
    hist_area.SetRange(0.0, 0.001, 1000);
    for (auto& f : m.face) {
        if (DistortionMetric::Area3D(f) > 0) {
            double areaUV = DistortionMetric::AreaUV(f);
            if (std::isfinite(areaUV)) {
                hist_area.Add(std::abs(areaUV), DistortionMetric::Area3D(f));
            }
        }
    }

    std::vector<std::vector<RasterizedParameterizationStats>> statsAtMipLevels
            = GetRasterizationStatsAtMipmapLevels(m, textureObject);

    std::vector<std::vector<GeometryImageStats>> geomStatsAtMipLevels
            = GetGeometryImageStatsAtMipmapLevels(m, textureObject);

    // occupancy
    long totalFragments = 0;
    long usedFragments = 0;
    for (auto& v : statsAtMipLevels) {
        totalFragments += (v[0].rw * v[0].rh);
        usedFragments += v[0].totalFragments - v[0].lostFragments;
    }
    double occupancy = usedFragments / (double) totalFragments;

    // write all the stats to a json object
    std::ofstream stats_file(filename + "_" + m.name + ".json");

    stats_file << "{" << std::endl;

    stats_file << JSONField("mesh"                , m.name)                        << "," << std::endl;
    stats_file << JSONField("vn"                  , m.FN())                        << "," << std::endl;
    stats_file << JSONField("fn"                  , m.VN())                        << "," << std::endl;
    stats_file << JSONField("connected_components", connected_components)          << "," << std::endl;
    stats_file << JSONField("boundary_loops"      , boundary_loops)                << "," << std::endl;
    stats_file << JSONField("genus"               , genus)                         << "," << std::endl;
    stats_file << JSONField("nonmanif_vert"       , nonmanif_vert)                 << "," << std::endl;
    stats_file << JSONField("nonmanif_edge"       , nonmanif_edge)                 << "," << std::endl;
    stats_file << JSONField("mapped_fraction"     , mapped_fraction)               << "," << std::endl;
    stats_file << JSONField("num_charts"          , num_charts)                    << "," << std::endl;
    stats_file << JSONField("num_null_charts"     , num_null_charts)               << "," << std::endl;
    stats_file << JSONField("occupancy"           , occupancy)                     << "," << std::endl;
    stats_file << JSONField("ntex"                , (int) statsAtMipLevels.size()) << "," << std::endl;

    for (std::size_t ntex = 0; ntex < statsAtMipLevels.size(); ntex++) {
        std::stringstream ss;
        ss << "texture_" << ntex;
        stats_file << JSONField(ss.str().c_str(), statsAtMipLevels[ntex]) << "," << std::endl;
    }

    for (std::size_t ntex = 0; ntex < geomStatsAtMipLevels.size(); ntex++) {
        std::stringstream ss;
        ss << "geometry_texture_" << ntex;
        stats_file << JSONField(ss.str().c_str(), geomStatsAtMipLevels[ntex]) << "," << std::endl;
    }

    stats_file << JSONField("qc_dist_avg"    , hist_angle.Avg()) << "," << std::endl;
    stats_file << JSONField("qc_dist_hist"   , hist_angle)       << "," << std::endl;
    stats_file << JSONField("qc_area_hist"   , hist_area)               << std::endl;

    stats_file << "}" << std::endl;

    stats_file.close();
}


void LogExecutionParameters(const Args& args, const ParameterizationStrategy& strategy)
{
    std::ofstream param_file("parameters.json");

    param_file << "{" << std::endl;
    param_file << JSONField("filename"           , args.filename)                << "," << std::endl;
    param_file << JSONField("regionCount"        , args.regionCount)             << "," << std::endl;
    param_file << JSONField("fixContextCapture"  , args.fixContextCapture)       << "," << std::endl;
    param_file << JSONField("optimizerIterations", strategy.optimizerIterations) << "," << std::endl;
    param_file << JSONField("scaffold"           , strategy.scaffold)                   << std::endl;
    param_file << "}" << std::endl;

    param_file.close();
}














