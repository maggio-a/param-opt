#ifndef MESH_UTILS_H
#define MESH_UTILS_H

#include "mesh.h"

#include <vector>

void MarkInitialSeamsAsCreases(Mesh& m, const SimpleTempData<TriMesh::FaceContainer, RegionID> &initialId);

void MarkShortestSeamToBorderAsNonFaux(Mesh& m, const face::Pos<Mesh::FaceType>& pos);

bool ComputePathToBoundary(Mesh& m, const face::Pos<Mesh::FaceType>& pos, std::vector<face::Pos<Mesh::FaceType>>& path);
void ClearCreasesAlongPath(Mesh& m, std::vector<PosF>& path);

#endif // MESH_UTILS_H

