#ifndef MESH_H
#define MESH_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>

#include <vector>
#include <memory>

#include <QImage>

#include "gl_utils.h"

using namespace vcg;

class MeshVertex;
class MeshFace;
//class MeshEdge;
struct MeshUsedTypes : public UsedTypes<Use<MeshVertex>::AsVertexType, Use<MeshFace>::AsFaceType/*, Use<MeshEdge>::AsEdgeType*/> {};

class MeshVertex : public Vertex<MeshUsedTypes, vertex::Coord3d, vertex::TexCoord2d, vertex::Normal3d,/* vertex::VEAdj,*/ vertex::VFAdj, vertex::Color4b, vertex::Qualityd, vertex::BitFlags> {};
class MeshFace : public Face<MeshUsedTypes, face::VertexRef, face::FFAdj, face::VFAdj, face::Mark, face::WedgeTexCoord2d, face::Normal3d, face::Color4b, face::Qualityf, face::BitFlags>
{
public:
    bool holeFilling = false;
    unsigned char isb = 0;

    void MarkSeamEdge(int i)
    {
        switch (i) {
        case 0:
            isb |= 1;
            break;
        case 1:
            isb |= 2;
            break;
        case 2:
            isb |= 4;
            break;
        default:
            assert(0 && "MarkSeamEdge");
        }
    }

    void UnmarkSeamEdge(int i)
    {
        switch (i) {
        case 0:
            isb &= ~1;
            break;
        case 1:
            isb &= ~2;
            break;
        case 2:
            isb &= ~4;
            break;
        default:
            assert(0 && "UnmarkSeamEdge");
        }
    }

    bool SeamEdge(int i)
    {
        assert(i >= 0 && i <= 2);
        return isb & (1 << i);
    }
};

//class MeshEdge : public Edge<MeshUsedTypes, edge::VertexRef, edge::VEAdj, edge::EEAdj, edge::BitFlags> {};
class Mesh : public tri::TriMesh<std::vector<MeshVertex>, std::vector<MeshFace>/*, std::vector<MeshEdge>*/> {};

using DistortionWedge = tri::Distortion<Mesh,true>;
using RegionID = std::size_t;

constexpr int INVALID_ID = 0xffffffff;

bool LoadMesh(Mesh &m, const char *fileName, TextureObjectHandle& textureObject, int &loadMask, std::string &modelName);
bool SaveMesh(Mesh &m, const char *fileName, TextureObjectHandle& textureObject, bool color = false);

// Builds a PMesh with face face topology and bounding box initialized, so that it can be
// passed to the poisson solver
void BuildMeshFromFacePointers(Mesh &m, const std::vector<std::vector<Mesh::FacePointer>* >& vFpVecp);

// Returns true if the mesh can be parameterized by the poisson solver
bool Parameterizable(Mesh &m);


#endif // MESH_H
