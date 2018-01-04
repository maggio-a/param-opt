#ifndef MESH_H
#define MESH_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>

#include <vector>
#include <memory>

#include <QImage>

#include "gl_util.h"

using namespace vcg;

class MeshVertex;
class MeshFace;
struct MeshUsedTypes : public UsedTypes<Use<MeshVertex>::AsVertexType, Use<MeshFace>::AsFaceType> {};

class MeshVertex : public Vertex<MeshUsedTypes, vertex::Coord3d, vertex::TexCoord2d, vertex::Normal3d, vertex::BitFlags> {};
class MeshFace : public Face<MeshUsedTypes, face::VertexRef, face::FFAdj, face::WedgeTexCoord2d, face::Color4b, face::Qualityf, face::BitFlags> {};
class Mesh : public tri::TriMesh<std::vector<MeshVertex>, std::vector<MeshFace>> {};

class PMeshVertex;
class PMeshFace;
struct PMeshUsedTypes : public UsedTypes<Use<PMeshVertex>::AsVertexType, Use<PMeshFace>::AsFaceType> {};

class PMeshVertex : public Vertex<PMeshUsedTypes, vertex::Coord3d, vertex::Normal3f, vertex::TexCoord2f, vertex::BitFlags> {};
class PMeshFace : public Face<PMeshUsedTypes, face::VertexRef, face::FFAdj, face::WedgeTexCoord2f, face::BitFlags> {};
class PMesh : public tri::TriMesh<std::vector<PMeshVertex>, std::vector<PMeshFace>> {};

using DistortionWedge = tri::Distortion<Mesh,true>;
using RegionID = std::size_t;

bool LoadMesh(Mesh &m, const char *fileName, TextureObjectHandle& textureObject, int &loadMask, std::string &modelName);
bool SaveMesh(Mesh &m, const char *fileName, TextureObjectHandle& textureObject);

// Builds a PMesh with face face topology and bounding box initialized, so that it can be
// passed to the poisson solver
template <typename MeshType>
static void BuildPMeshFromFacePointers(PMesh &pm, const std::vector<std::vector<typename MeshType::FacePointer>* >& vFpVecp)
{
    pm.Clear();

    auto f = [&pm](typename MeshType::FacePointer fptr) {
        tri::Allocator<PMesh>::AddFace(pm, fptr->P(0), fptr->P(1), fptr->P(2));
    };

    for (auto fpVecp : vFpVecp) std::for_each(fpVecp->begin(), fpVecp->end(), f);

    tri::Clean<PMesh>::RemoveDuplicateVertex(pm);
    tri::Allocator<PMesh>::CompactEveryVector(pm);

    tri::UpdateTopology<PMesh>::FaceFace(pm);
    tri::UpdateBounding<PMesh>::Box(pm);
}

// Returns true if the mesh can be parameterized by the poisson solver
template <typename MeshType>
static bool Parameterizable(MeshType &m)
{
    if (tri::Clean<MeshType>::CountNonManifoldEdgeFF(m) > 0) {
        return false;
    }

    if (tri::Clean<MeshType>::CountNonManifoldVertexFF(m) > 0) {
        return false;
    }

    if (tri::Clean<MeshType>::IsWaterTight(m)) {
        return false;
    }

    if (tri::Clean<MeshType>::MeshGenus(m) > 0) {
        return false;
    }

    return true;
}

#endif // MESH_H
