#ifndef TEXTURE_RENDERING_H
#define TEXTURE_RENDERING_H

#include "gl_utils.h"
#include "mesh_graph.h"

#include <vector>

class Mesh;
class MeshFace;
class GLFWwindow;

enum InterpolationMode {
    Nearest, Linear, Cubic
};

struct TextureSize {
    int w;
    int h;
};

struct RasterizedParameterizationStats {
    int rw; // raster width
    int rh; // raster height
    int totalFragments;
    int totalFragments_bilinear;
    int overwrittenFragments; // number of fragments that were written more than once
    int lostFragments; // number of fragments lost due to overwrites: if fragment f has fw>1 writes, than lostFragmens += (fw-1)
    int boundaryFragments;
};

std::vector<TextureSize> ComputeSizes(int nTex, TextureObjectHandle inputTexture);

TextureObjectHandle RenderTexture(Mesh &m, TextureObjectHandle textureObject, bool filter, InterpolationMode imode, GLFWwindow *parentWindow);

RasterizedParameterizationStats GetRasterizationStats(ChartHandle chart, int width, int height);

std::vector<RasterizedParameterizationStats> GetRasterizationStats(Mesh& m, TextureObjectHandle textureObject);

#endif // TEXTURE_RENDERING_H

