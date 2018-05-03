#ifndef MESH_ATTRIBUTE_H
#define MESH_ATTRIBUTE_H

#include "mesh.h"

struct TexCoordStorage {
    vcg::TexCoord2d tc[3];
};

struct CoordStorage {
    vcg::Point3d P[3];
};

struct ParameterizationScaleInfo {
    double surfaceArea;
    double parameterArea;
    double scale;
    int numNonZero;

    ParameterizationScaleInfo() : surfaceArea{0}, parameterArea{0}, scale{0}, numNonZero{0} {}
};

inline Mesh::PerMeshAttributeHandle<ParameterizationScaleInfo> GetParameterizationScaleInfoAttribute(Mesh& m)
{
    return tri::Allocator<Mesh>::GetPerMeshAttribute<ParameterizationScaleInfo>(m, "MeshAttribute_ParameterizationScaleInfo");
}

inline Mesh::PerFaceAttributeHandle<TexCoordStorage> GetWedgeTexCoordStorageAttribute(Mesh& m)
{
    return tri::Allocator<Mesh>::GetPerFaceAttribute<TexCoordStorage>(m, "WedgeTexCoordStorage");
}

inline Mesh::PerFaceAttributeHandle<CoordStorage> GetTargetShapeAttribute(Mesh& shell)
{
    return tri::Allocator<Mesh>::GetPerFaceAttribute<CoordStorage>(shell, "FaceAttribute_TargetShape");
}

/* The index of the initial mesh face that corresponds to each shell face
 * Note that each shell face is either mapped to a mesh face, or it is hole-filling */
inline Mesh::PerFaceAttributeHandle<int> GetFaceIndexAttribute(Mesh& shell)
{
    return tri::Allocator<Mesh>::GetPerFaceAttribute<int>(shell, "FaceAttribute_FaceIndex");
}


#endif // MESH_ATTRIBUTE_H

