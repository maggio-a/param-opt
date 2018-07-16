#ifndef MESH_H
#define MESH_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>

#include <vector>
#include <memory>

#include <QImage>

#include "gl_utils.h"

using namespace vcg;

/* per face extra flags, defined according to the vcg style */
template <class T>
class FaceQualifier : public T {

public:

    FaceQualifier() : _qualifier(0) {}

    void SetMesh()        { _qualifier = MESH; }
    void SetHoleFilling() { _qualifier = HOLE_FILLING; }
    void SetScaffold()    { _qualifier = SCAFFOLD; }

    bool IsMesh()        const { return _qualifier == MESH; }
    bool IsHoleFilling() const { return _qualifier == HOLE_FILLING; }
    bool IsScaffold()    const { return _qualifier == SCAFFOLD; }

    template <class RType>
    void ImportData(const RType& rhs) {
        rhs.HasQualifier() && "RHS NOT COMPATIBLE";
        _qualifier = rhs._qualifier;
        T::ImportData(rhs);
    }
    void Alloc(const int & ns){ T::Alloc(ns); }
    void Dealloc(){ T::Dealloc(); }
    static bool HasQualifier() { return true; }
    static void Name(std::vector<std::string> & name) {
        name.push_back(std::string("FaceQualifier"));
        T::Name(name);
    }

private:
    static constexpr unsigned char MESH = 1;
    static constexpr unsigned char HOLE_FILLING = 2;
    static constexpr unsigned char SCAFFOLD = 3;

    unsigned char _qualifier;
};

class MeshVertex;
class MeshFace;
//class MeshEdge;
struct MeshUsedTypes : public UsedTypes<Use<MeshVertex>::AsVertexType, Use<MeshFace>::AsFaceType/*, Use<MeshEdge>::AsEdgeType*/> {};

class MeshVertex : public Vertex<MeshUsedTypes, vertex::Coord3d, vertex::TexCoord2d, vertex::Normal3d,/* vertex::VEAdj,*/ vertex::VFAdj, vertex::Color4b, vertex::Qualityd, vertex::BitFlags> {};
class MeshFace : public Face<MeshUsedTypes, FaceQualifier, face::VertexRef, face::FFAdj, face::VFAdj, face::Mark, face::WedgeTexCoord2d, face::Normal3d, face::Color4b, face::Qualityf, face::BitFlags>
{
};

//class MeshEdge : public Edge<MeshUsedTypes, edge::VertexRef, edge::VEAdj, edge::EEAdj, edge::BitFlags> {};
class Mesh : public tri::TriMesh<std::vector<MeshVertex>, std::vector<MeshFace>/*, std::vector<MeshEdge>*/>{
public:
    std::string name{"mesh"};
};


class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>::AsVertexType, Use<MyEdge>::AsEdgeType, Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3d, vertex::Normal3d, vertex::VEAdj, vertex::VFAdj,vertex::BitFlags  >{};
class MyEdge    : public Edge<   MyUsedTypes, edge::VertexRef, edge::VEAdj,     edge::EEAdj, edge::BitFlags> {};
class MyFace    : public Face  < MyUsedTypes, face::VertexRef, face::VFAdj, face::FFAdj, face::Mark, face::Color4b, face::BitFlags > {};
class MyMesh : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyEdge>, std::vector<MyFace> >{};


//using DistortionWedge = tri::Distortion<Mesh,true>;
using RegionID = std::size_t;

constexpr RegionID INVALID_ID = 0xffffffff;


bool LoadMesh(Mesh &m, const char *fileName, TextureObjectHandle& textureObject, int &loadMask);
bool SaveMesh(Mesh &m, const char *fileName, TextureObjectHandle& textureObject, bool color = false);

/* Builds a mesh from a given vector of face pointers. The order of the faces
 * is guaranteed to be preserved in the face container of the mesh. */
void MeshFromFacePointers(const std::vector<Mesh::FacePointer>& vfp, Mesh& out);

// Builds a PMesh with face face topology and bounding box initialized, so that it can be
// passed to the poisson solver
void BuildMeshFromFacePointers(Mesh &m, const std::vector<std::vector<Mesh::FacePointer>* >& vFpVecp);

// Returns true if the mesh can be parameterized by the poisson solver
bool Parameterizable(Mesh &m);


#endif // MESH_H
