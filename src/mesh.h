#ifndef MESH_H
#define MESH_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>

#include <vector>
#include <memory>

#include <QImage>

using namespace vcg;

class MeshVertex;
class MeshFace;
struct MeshUsedTypes : public UsedTypes<Use<MeshVertex>::AsVertexType, Use<MeshFace>::AsFaceType> {};

class MeshVertex : public Vertex<MeshUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags> {};
class MeshFace : public Face<MeshUsedTypes, face::VertexRef, face::FFAdj, face::WedgeTexCoord2f, face::BitFlags> {};
class Mesh : public tri::TriMesh<std::vector<MeshVertex>, std::vector<MeshFace>> {};

class PMeshVertex;
class PMeshFace;
struct PMeshUsedTypes : public UsedTypes<Use<PMeshVertex>::AsVertexType, Use<PMeshFace>::AsFaceType> {};

class PMeshVertex : public Vertex<PMeshUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::TexCoord2f, vertex::BitFlags> {};
class PMeshFace : public Face<PMeshUsedTypes, face::VertexRef, face::FFAdj, face::WedgeTexCoord2f, face::BitFlags> {};
class PMesh : public tri::TriMesh<std::vector<PMeshVertex>, std::vector<PMeshFace>> {};

using DistortionWedge = tri::Distortion<Mesh,true>;

using RegionID = std::size_t;

bool LoadMesh(Mesh &m, const char *fileName, std::vector<std::shared_ptr<QImage>> &imgVec, int &loadMask, std::string &modelName);

#endif // MESH_H
