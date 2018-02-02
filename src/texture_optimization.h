#ifndef TEXTURE_OPTIMIZATION_H
#define TEXTURE_OPTIMIZATION_H

#include "mesh_graph.h"

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

/// TODO REFACTOR (these enums should be each in their respective headers)
enum DirectParameterizer {
    DCP, FixedBorderBijective
};

enum TexCoordOptimizer {
    AreaPreserving, SymmetricDirichletOpt, MIPS
};

enum DescentType {
    Gradient, LimitedMemoryBFGS, ScalableLocallyInjectiveMappings
};

struct ParameterizationStrategy {
    DirectParameterizer directParameterizer = DCP;
    TexCoordOptimizer optimizer = SymmetricDirichletOpt;
    ParameterizationGeometry geometry = Model;
    DescentType descent = Gradient;
    int optimizerIterations = 0;
};

void ReparameterizeZeroAreaRegions(Mesh &m, std::shared_ptr<MeshGraph> graph);
void PreprocessMesh(Mesh& m);

bool ChartParameterizationHasOverlaps(Mesh& m, GraphManager::ChartHandle chart);

bool ParameterizeMesh(Mesh& m, ParameterizationStrategy strategy);

bool ParameterizeChart(Mesh &m, GraphManager::ChartHandle ch, ParameterizationStrategy strategy);

/*
 * Parameterize the mesh graph. Each region is a group of connected mesh faces, and it is assumed to be homeomorphic to a disk.
 * The parameterization of each region is performed according to the parameterization strategy passed as parameter.
 * If failsafe is true, then each region is tested for overlaps in the parameterization: in case the parameterization contains
 * overlaps, the region is split in its original components and each face is assigned its original texture coordinates. This
 * ensures that no overlaps are introduced by this procedure, potentially reducing the whole procedure to a no-op if necessary
 * (that is, if every region parameterization contains overlaps).
 * The threshold parameter is the fraction of overlapping fragments - in the rasterized parameterization - above which the chart
 * is split.
 *
 * After each region is parameterized the procedure packs the texture atlas.
 *
 * Returns the number of charts that could not be parameterized.
 */
/// TODO update distortion info if needed (this should also be done through the graph manager)
int ParameterizeGraph(GraphManager& gm,
                      ParameterizationStrategy strategy,
                      bool failsafe,
                      double threshold = 0);


void ReduceTextureFragmentation_NoPacking(GraphManager &gm, std::size_t minRegionSize);
void ReduceTextureFragmentation(Mesh &m, std::shared_ptr<MeshGraph> graph, std::size_t minRegionSize);

#endif // TEXTURE_OPTIMIZATION_H
