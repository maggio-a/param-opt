#ifndef UV_H
#define UV_H

#include "mesh.h"
#include <vcg/space/texcoord2.h>
#include <vcg/space/box2.h>

/* Save a copy of the original texture coordinates (this will be used to render the new texture) */
void StoreWedgeTexCoordAsAttribute(Mesh &m);

vcg::Box2d UVBox(const Mesh& m);

/* Computes and stores as attribute the scale from 3D to UV space of the existing parameterization */
double ComputeParameterizationScaleInfo(Mesh& m);

/* Computes per face connected component ids. Two faces with the same id belong
 * to the same connected component. Uses FaceFaceFromTexCoord topology.
 * */
std::size_t ComputePerFaceConnectedComponentIdAttribute(Mesh &m);

#endif // UV_H
