#ifndef LOGGING_H
#define LOGGING_H

#include <vcg/math/histogram.h>


#include <memory>

struct ParameterizationStrategy;
struct MeshGraph;
struct RasterizedParameterizationStats;

void LogStrategy(const ParameterizationStrategy& strategy, double tol);
void PercentilePlot(vcg::Distribution<double>& d);
void LogDistortionStats(std::shared_ptr<MeshGraph> graph);
void LogParameterizationStats(std::shared_ptr<MeshGraph> graph, const std::vector<RasterizedParameterizationStats>& statsVec);

#endif // LOGGING_H

