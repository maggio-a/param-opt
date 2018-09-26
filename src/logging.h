#ifndef LOGGING_H
#define LOGGING_H

#include <vcg/math/histogram.h>

#include <string>
#include <memory>

#include "gl_utils.h"

struct ParameterizationStrategy;
struct MeshGraph;
struct RasterizedParameterizationStats;
struct Args;

void LogStrategy(const ParameterizationStrategy& strategy, double tol);
void PercentilePlot(vcg::Distribution<double>& d);
void LogDistortionStats(std::shared_ptr<MeshGraph> graph);
void LogParameterizationStats(std::shared_ptr<MeshGraph> graph, const std::vector<RasterizedParameterizationStats>& statsVec);

void LogAggregateStats(const std::string& str, std::shared_ptr<MeshGraph> graph, TextureObjectHandle textureObject);

void LogExecutionParameters(const Args& args, const ParameterizationStrategy& strategy);

#endif // LOGGING_H

