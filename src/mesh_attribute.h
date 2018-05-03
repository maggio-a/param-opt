#ifndef MESH_ATTRIBUTE_H
#define MESH_ATTRIBUTE_H

#include "mesh.h"

/*
 * The input mesh has the following attributes defined
 *
 * PER MESH
 *
 *   ParameterizationScaleInfo (GetParameterizationScaleInfoAttribute)
 *
 * PER FACE
 *
 *   - WedgeTexCoordStorage (GetWedgeTexCoordStorageAttribute)
 *       The per-wedge texture coordinates of the input mesh
 *
 *   - ConnectedComponentID (GetConnectedComponentIDAttribute)
 *       The id of the segment to which the face belongs. The idea is that a segment is a set of connected faces, that are parameterized
 *       together and will end up in the same chart in the texture atlas
 *
 *   - InitialConnectedComponentID (GetInitialConnectedComponentIDAttribute)
 *       The id of the chart a face belonged in the input parameterization.
 *
 *
 *
 * A shell mesh has the following attributes defined
 *
 * PER FACE
 *
 *   - FaceAttribute_TargetShape (GetTargetShapeAttribute)
 *       The shape of the triangle towards which distortion must be minimized when the shell is optimized.
 *
 *   - FaceAttribute_FaceIndex (GetFaceIndexAttribute)
 *       The index of the input mesh face that corresponds to the shell face. For faces added to fill holes,
 *       the index is -1
 * */

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

inline Mesh::PerFaceAttributeHandle<RegionID> GetConnectedComponentIDAttribute(Mesh& m)
{
    return tri::Allocator<Mesh>::GetPerFaceAttribute<RegionID>(m, "ConnectedComponentID");
}

inline bool HasConnectedComponentIDAttribute(Mesh& m)
{
    return tri::Allocator<Mesh>::IsValidHandle<RegionID>(m, tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "ConnectedComponentID"));
}

inline Mesh::PerFaceAttributeHandle<RegionID> GetInitialConnectedComponentIDAttribute(Mesh& m)
{
    return tri::Allocator<Mesh>::GetPerFaceAttribute<RegionID>(m, "InitialConnectedComponentID");
}

inline bool HasInitialConnectedComponentIDAttribute(Mesh& m)
{
    return tri::Allocator<Mesh>::IsValidHandle<RegionID>(m, tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "InitialConnectedComponentID"));
}

inline Mesh::PerFaceAttributeHandle<CoordStorage> GetTargetShapeAttribute(Mesh& shell)
{
    return tri::Allocator<Mesh>::GetPerFaceAttribute<CoordStorage>(shell, "FaceAttribute_TargetShape");
}

inline bool HasTargetShapeAttribute(Mesh& shell)
{
    return tri::Allocator<Mesh>::IsValidHandle<RegionID>(shell, tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(shell, "FaceAttribute_TargetShape"));
}

inline Mesh::PerFaceAttributeHandle<int> GetFaceIndexAttribute(Mesh& shell)
{
    return tri::Allocator<Mesh>::GetPerFaceAttribute<int>(shell, "FaceAttribute_FaceIndex");
}

#endif // MESH_ATTRIBUTE_H

