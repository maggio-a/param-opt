#ifndef MESH_ATTRIBUTE_H
#define MESH_ATTRIBUTE_H

#include "mesh.h"

/*
 * The input mesh has the following attributes defined
 *
 * PER MESH
 *
 *   - ParameterizationScaleInfo (GetParameterizationScaleInfoAttribute)
 *       The scale factor from 3D to UV. Only parameterized faces are taken into
 *       account
 *
 * PER FACE
 *
 *   - WedgeTexCoordStorage (GetWedgeTexCoordStorageAttribute)
 *       The per-wedge texture coordinates of the input mesh
 *
 *   - ConnectedComponentID (GetConnectedComponentIDAttribute)
 *       The id of the segment to which the face belongs. The idea is that a
 *       segment is a set of connected faces, that are parameterized together
 *       and will end up in the same chart in the texture atlas
 *
 *   - InitialConnectedComponentID (GetInitialConnectedComponentIDAttribute)
 *       The id of the chart a face belonged in the input parameterization.
 *
 *
 *
 * A shell mesh has the following attributes defined
 *
 * PER MESH
 *
 *   - BoundaryInfo (GetBoundaryInfoAttribute)
 *       The BoundaryInfo structure holds information about the boundary of the
 *       shell in 3D, this is computed when the shell is created, and then is
 *       never updated after the shell shape is projected to parameter space. In
 *       doing so the boundary information is 'fixed' and not influenced by the
 *       shell parameterization/syncing operations.
 *
 * PER FACE
 *
 *   - FaceAttribute_TargetShape (GetTargetShapeAttribute)
 *       The shape of the triangle towards which distortion must be minimized
 *       when the shell is optimized.
 *
 *   - FaceAttribute_FaceIndex (GetFaceIndexAttribute)
 *       The index of the input mesh face that corresponds to the shell face.
 *       For faces added to fill holes, the index is -1
 *
 *   - FaceAttribute_HoleFilling (GetHoleFillingAttribute)
 *       A boolean indicating whether the shell face is added to fill holes or
 *       is a genuine face.
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

struct BoundaryInfo {
    std::vector<double> vBoundaryLength;
    std::vector<std::size_t> vBoundarySize;
    std::vector<std::vector<std::size_t>> vBoundaryFaces;
    std::vector<std::vector<int>> vVi; // The face boundary vertex indices

    std::size_t N();
    std::size_t LongestBoundary();
    void Clear();
};

inline std::size_t BoundaryInfo::N()
{
    assert(vBoundaryLength.size() == vBoundarySize.size() && vBoundaryLength.size() == vBoundaryFaces.size());
    return vBoundaryLength.size();
}

inline std::size_t BoundaryInfo::LongestBoundary()
{
    assert(N() > 0);
    return std::distance(vBoundaryLength.begin(),
                         std::max_element(vBoundaryLength.begin(), vBoundaryLength.end()));
}

inline void BoundaryInfo::Clear()
 {
    vBoundaryLength.clear();
    vBoundarySize.clear();
    vBoundaryFaces.clear();
    vVi.clear();
}

inline Mesh::PerMeshAttributeHandle<ParameterizationScaleInfo> GetParameterizationScaleInfoAttribute(Mesh& m)
{
    return tri::Allocator<Mesh>::GetPerMeshAttribute<ParameterizationScaleInfo>(m, "MeshAttribute_ParameterizationScaleInfo");
}

inline bool HasParameterizationScaleInfoAttribute(Mesh& m)
{
    return tri::Allocator<Mesh>::IsValidHandle<ParameterizationScaleInfo>(
                m, tri::Allocator<Mesh>::FindPerMeshAttribute<ParameterizationScaleInfo>(m, "MeshAttribute_ParameterizationScaleInfo"));
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
    return tri::Allocator<Mesh>::IsValidHandle<RegionID>(
                m, tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "ConnectedComponentID"));
}

inline Mesh::PerFaceAttributeHandle<RegionID> GetInitialConnectedComponentIDAttribute(Mesh& m)
{
    return tri::Allocator<Mesh>::GetPerFaceAttribute<RegionID>(m, "InitialConnectedComponentID");
}

inline bool HasInitialConnectedComponentIDAttribute(Mesh& m)
{
    return tri::Allocator<Mesh>::IsValidHandle<RegionID>(
                m, tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "InitialConnectedComponentID"));
}

inline Mesh::PerMeshAttributeHandle<BoundaryInfo> GetBoundaryInfoAttribute(Mesh& shell)
{
    return tri::Allocator<Mesh>::GetPerMeshAttribute<BoundaryInfo>(shell, "MeshAttribute_BoundaryInfo");
}

inline bool HasBoundaryInfoAttribute(Mesh& shell)
{
    return tri::Allocator<Mesh>::IsValidHandle<BoundaryInfo>(
                shell, tri::Allocator<Mesh>::FindPerMeshAttribute<BoundaryInfo>(shell, "MeshAttribute_BoundaryInfo"));
}

inline Mesh::PerFaceAttributeHandle<CoordStorage> GetTargetShapeAttribute(Mesh& shell)
{
    return tri::Allocator<Mesh>::GetPerFaceAttribute<CoordStorage>(shell, "FaceAttribute_TargetShape");
}

inline bool HasTargetShapeAttribute(Mesh& shell)
{
    return tri::Allocator<Mesh>::IsValidHandle<CoordStorage>(
                shell, tri::Allocator<Mesh>::FindPerFaceAttribute<CoordStorage>(shell, "FaceAttribute_TargetShape"));
}

inline Mesh::PerFaceAttributeHandle<int> GetFaceIndexAttribute(Mesh& shell)
{
    return tri::Allocator<Mesh>::GetPerFaceAttribute<int>(shell, "FaceAttribute_FaceIndex");
}

inline bool HasFaceIndexAttribute(Mesh& shell)
{
    return tri::Allocator<Mesh>::IsValidHandle<int>(
                shell, tri::Allocator<Mesh>::FindPerFaceAttribute<int>(shell, "FaceAttribute_FaceIndex"));
}

#endif // MESH_ATTRIBUTE_H

