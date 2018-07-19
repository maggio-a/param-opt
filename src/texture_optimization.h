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
void PreprocessMesh(Mesh& m, GraphHandle graph);

bool ChartParameterizationHasOverlaps(Mesh& m, GraphManager::ChartHandle chart);

/* Parameterize the mesh graph. The parameterization of each region is performed
 * according to the ParameterizationStrategy passed.
 * The injectivityTolerance parameter is the fraction of overlapping fragments
 * - in the rasterized parameterization - above which backtracking is triggered.
 * Returns the number of charts that could not be parameterized. */
/// TODO update distortion info if needed (this should also be done through the graph manager)
int ParameterizeGraph(GraphManager& gm, ParameterizationStrategy strategy, double injectivityTolerance);

void RecoverFromSplit(std::vector<ChartHandle>& split, GraphManager& gm, std::vector<ChartHandle>& chartVec, bool binarySplit);

/* Pack the texture atlas encoded in the graph. Assumes the segmentation
 * correctly reflects the thexture coordinates */
void Pack(GraphHandle graph);


void RecomputeSegmentation(GraphManager &gm, std::size_t regionCount, std::size_t minRegionSize);

#endif // TEXTURE_OPTIMIZATION_H
