#ifndef MESH_H
#define MESH_H

#include <vcg/complex/complex.h>

#include <vector>

using namespace vcg;

class MeshVertex;
class MeshFace;
struct MeshUsedTypes : public UsedTypes<Use<MeshVertex>::AsVertexType, Use<MeshFace>::AsFaceType> {};

class MeshVertex : public Vertex<MeshUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags> {};
class MeshFace : public Face<MeshUsedTypes, face::VertexRef, face::FFAdj, face::WedgeTexCoord2f, face::BitFlags> {};
class Mesh : public tri::TriMesh<std::vector<MeshVertex>, std::vector<MeshFace>> {};

struct WedgeTexCoordStorage {
    TexCoord2f tc[3];
};

#endif // MESH_H

