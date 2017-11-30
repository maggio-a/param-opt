#ifndef VERTEX_POSITION_H
#define VERTEX_POSITION_H

#include <cassert>

template <typename MeshType>
struct DefaultVertexPosition
{
    using FacePointer = typename MeshType::FacePointer;
    using CoordType = typename MeshType::CoordType;

    CoordType operator()(FacePointer fp, int i) { return fp->V(i)->P(); }

    CoordType V0(FacePointer fp, int i = 0) { return (*this)(fp, i); }
    CoordType V1(FacePointer fp, int i = 0) { return (*this)(fp, (i+1)%3); }
    CoordType V2(FacePointer fp, int i = 0) { return (*this)(fp, (i+2)%3); }
};

template <typename MeshType>
struct WedgeTexCoordVertexPosition
{
    using FacePointer = typename MeshType::FacePointer;
    using CoordType = typename MeshType::CoordType;

    CoordType operator()(FacePointer fp, int i) { assert(i>=0 && i<3); return CoordType{fp->WT(i).U(), fp->WT(i).V(), 0}; }

    CoordType V0(FacePointer fp, int i = 0) { return (*this)(fp, i); }
    CoordType V1(FacePointer fp, int i = 0) { return (*this)(fp, (i+1)%3); }
    CoordType V2(FacePointer fp, int i = 0) { return (*this)(fp, (i+2)%3); }
};

#endif // VERTEX_POSITION_H

