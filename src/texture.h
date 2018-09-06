#ifndef TEXTURE_H
#define TEXTURE_H

#include "mesh.h"
#include "mesh_graph.h"
#include "gl_utils.h"

#include <unordered_map>

void GenerateDistortionTextures(Mesh& m, TextureObjectHandle textureObject);

void GeneratePackingQualityTexture(GraphHandle graph, TextureObjectHandle textureObject, std::unordered_map<RegionID, float> &qualityMap);

void EvaluateGutterResistance(Mesh& m, TextureObjectHandle textureObject);

#endif // TEXTURE_H

