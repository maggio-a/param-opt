#ifndef TEXTURE_OPTIMIZATION_H
#define TEXTURE_OPTIMIZATION_H

#include "mesh_graph.h"
#include "iterative.h"
#include "parameterization.h"

#include <vcg/space/point2.h>

#include <utility>

class Mesh;

struct Point2iHasher {
    std::size_t operator()(const Point2i& p) const noexcept
    {
        std::size_t seed = 0;
        seed ^= std::hash<int>()(p[0]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<int>()(p[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

enum TexCoordOptimizer {
    AreaPreserving, SymmetricDirichletOpt, MIPS
};

void ReparameterizeZeroAreaRegions(Mesh &m, std::shared_ptr<MeshGraph> graph);
void PreprocessMesh(Mesh& m);

bool ChartParameterizationHasOverlaps(Mesh& m, GraphManager::ChartHandle chart);

/*
 * Parameterize the mesh graph. The parameterization of each region is performed
 * according to the ParameterizationStrategy passed.
 * The injectivityTolerance parameter is the fraction of overlapping fragments
 * in the rasterized parameterization above which the chart is split into its
 * original components.
 * If retry is true, the original components are merged into two sub-charts and
 * parameterized again, until valid parameterizations are produced or the split
 * can no longer be performed, in which case the original texture coordinates
 * are restored.
 * After each region is parameterized the procedure packs the texture atlas.
 * Returns the number of charts that could not be parameterized. */
/// TODO update distortion info if needed (this should also be done through the graph manager)
int ParameterizeGraph(GraphManager& gm, ParameterizationStrategy strategy, double injectivityTolerance, bool retry);


void RecomputeSegmentation(GraphManager &gm, std::size_t regionCount, std::size_t minRegionSize);

#endif // TEXTURE_OPTIMIZATION_H
