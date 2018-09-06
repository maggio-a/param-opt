/*
 * References for the bicubic interpolated texture lookup:
 *  - GPU gems 2 ch 20 (Sigg and Hadwiger 2005)
 *  - Efficient GPU-Based Texture Interpolation using Uniform B-Splines  (Ruijters et al. 2009)
 * */

#ifndef TEXTURE_RENDERING_H
#define TEXTURE_RENDERING_H

#include "gl_utils.h"

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

std::vector<TextureSize> ComputeSizes(int nTex, TextureObjectHandle inputTexture);

TextureObjectHandle RenderTexture(Mesh &m, TextureObjectHandle textureObject, bool filter, InterpolationMode imode, GLFWwindow *parentWindow);


#endif // TEXTURE_RENDERING_H

