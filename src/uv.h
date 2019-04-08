#ifndef UV_H
#define UV_H

#include <vcg/space/box2.h>

#include "gl_utils.h"

class Mesh;
class MeshFace;

/*
 * Note on the management of texture coordinates.
 * To simplify things, texture coordinates are scaled relative to the texture
 * image size. This has two advantages: first, it naturally models the difference
 * in resolution of charts that have comparable areas but belong texture files
 * of different resolution. Second, we do not have to scale according to the ratio
 * of the texture dimensions when dealing with coordinates for rectangular textures,
 * which was necessary when computing for example the angle distortion of the
 * uv-to-3d mapping
 */

void ScaleTextureCoordinatesToImage(Mesh& m, TextureObjectHandle textureObject);

void ScaleTextureCoordinatesToParameterArea(Mesh& m, TextureObjectHandle textureObject);

/* Save a copy of the original texture coordinates (this will be used to render
 * the new texture) */
void StoreWedgeTexCoordAsAttribute(Mesh& m);

vcg::Box2d UVBox(const Mesh& m);
vcg::Box2d UVBoxVertex(const Mesh& m);

/* Computes and stores as attribute the scale factor from 3D to UV space of the
 * existing parameterization */
double ComputeParameterizationScaleInfo(Mesh& m);

/* Computes per face connected component ids. Two faces with the same id belong
 * to the same connected component. Uses FaceFaceFromTexCoord topology. */
std::size_t ComputePerFaceConnectedComponentIdAttribute(Mesh &m);

/* Marks the texture seams of the mesh (computed using FaceFaceFromTexCoord
 * topology as faux edges */
void MarkSeamsAsFaux(Mesh& m);

bool CheckLocalInjectivity(Mesh& m);

bool IsBorderUV(MeshFace *f, int e);

#endif // UV_H
