#ifndef VERTEX_POSITION_H
#define VERTEX_POSITION_H

#include <cassert>

#include "uv.h"

template <typename MeshType>
struct DefaultVertexPosition
{
    using FacePointer = typename MeshType::FacePointer;
    using ConstFacePointer = typename MeshType::ConstFacePointer;
    using CoordType = typename MeshType::CoordType;

    CoordType operator()(ConstFacePointer fp, int i) { assert(i>=0 && i<3); return fp->cV(i)->P(); }

    CoordType V0(ConstFacePointer fp, int i = 0) { return (*this)(fp, i); }
    CoordType V1(ConstFacePointer fp, int i = 0) { return (*this)(fp, (i+1)%3); }
    CoordType V2(ConstFacePointer fp, int i = 0) { return (*this)(fp, (i+2)%3); }
};

template <typename MeshType>
struct WedgeTexCoordAttributePosition
{
    using FacePointer = typename MeshType::FacePointer;
    using ConstFacePointer = typename MeshType::ConstFacePointer;
    using CoordType = typename MeshType::CoordType;

    typename MeshType::template PerFaceAttributeHandle<TexCoordStorage> attr;

    WedgeTexCoordAttributePosition(MeshType& m, const char *attributeName)
    {
        attr = tri::Allocator<MeshType>::template FindPerFaceAttribute<TexCoordStorage>(m, attributeName);
        assert(tri::Allocator<MeshType>::template IsValidHandle<TexCoordStorage>(m, attr));
    }

    CoordType operator()(ConstFacePointer fp, int i) { assert(i>=0 && i<3); return CoordType(attr[fp].tc[i].U(), attr[fp].tc[i].V(), 0); }

    CoordType V0(ConstFacePointer fp, int i = 0) { return (*this)(fp, i); }
    CoordType V1(ConstFacePointer fp, int i = 0) { return (*this)(fp, (i+1)%3); }
    CoordType V2(ConstFacePointer fp, int i = 0) { return (*this)(fp, (i+2)%3); }
};

#endif // VERTEX_POSITION_H

