#ifndef MESH_UTILS_H
#define MESH_UTILS_H

#include "mesh.h"

#include <vector>

using PosF = face::Pos<MeshFace>;

std::vector<PosF> GetCreasePosFan(PosF& startPos);
void MarkInitialSeamsAsCreases(Mesh& m, SimpleTempData<Mesh::FaceContainer, RegionID> &initialId);
void ClearCreaseLoops(Mesh& m);

void MarkShortestSeamToBorderAsNonFaux(Mesh& m, const PosF& pos);

double ComputeDistanceFromBorderOnSeams(Mesh& m);
bool ComputePathToBoundary(Mesh& m, const PosF& pos, std::vector<PosF>& path);
void ClearCreasesAlongPath(Mesh& m, std::vector<PosF>& path);

#endif // MESH_UTILS_H

