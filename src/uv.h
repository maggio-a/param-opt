#ifndef UV_H
#define UV_H

#include "mesh.h"
#include <vcg/space/texcoord2.h>
#include <vcg/space/box2.h>


struct TexCoordStorage {
    vcg::TexCoord2d tc[3];
};

struct CoordStorage {
    vcg::Point3d P[3];
};

/* Save a copy of the original texture coordinates (this will be used to render the new texture) */
void StoreWedgeTexCoordAsAttribute(Mesh &m);

vcg::Box2d UVBox(const Mesh& m);

/* Computes and stores as attribute the scale from 3D to UV space of the existing parameterization */
double ComputeParameterizationScaleInfo(Mesh& m);

/* Computes per face connected component ids. Two faces with the same id belong
 * to the same connected component. Uses FaceFaceFromTexCoord topology.
 * */
template <class MeshType>
std::size_t ComputePerFaceConnectedComponentIdAttribute(MeshType &m)
{
    typename MeshType::template PerFaceAttributeHandle<RegionID> CCIDh
            = tri::Allocator<MeshType>::template GetPerFaceAttribute<RegionID>(m, "ConnectedComponentID");

    tri::UpdateFlags<MeshType>::FaceClearV(m);
    tri::UpdateTopology<MeshType>::FaceFaceFromTexCoord(m);

    std::stack<typename MeshType::FacePointer> s;
    size_t regionCounter = 0;

    for (auto &f : m.face) {
        if (!f.IsV()) {
            f.SetV();
            s.push(&f);
            size_t id = regionCounter++;
            while (!s.empty()) {
                auto fp = s.top();
                s.pop();
                CCIDh[fp] = id;
                for (int i = 0; i < fp->VN(); ++i) {
                    if (!face::IsBorder(*fp, i)) {
                        auto adj = fp->FFp(i);
                        if (!adj->IsV()) {
                            adj->SetV();
                            s.push(adj);
                        }
                    }
                }

            }
        }
    }
    return regionCounter;
}

template<class ScalarType, class MeshType>
std::size_t ConvertTextureBoundaryToOutline2Vec(MeshType &m, std::vector<std::vector<Point2<ScalarType>>> &outline2Vec)
{
    typedef typename MeshType::FaceType FaceType;

    tri::RequirePerFaceWedgeTexCoord(m);
    std::vector<Point2<ScalarType>> outline;

    //tri::Allocator<MeshType>::CompactVertexVector(m);
    //tri::Allocator<MeshType>::CompactFaceVector(m);
    tri::UpdateFlags<MeshType>::FaceClearV(m);
    tri::UpdateFlags<MeshType>::VertexClearV(m);
    tri::UpdateTopology<MeshType>::FaceFace(m);

    for(auto &f : m.face) {
        for (int j=0; j < f.VN(); j++) {
            if (!f.IsV() && face::IsBorder(f, j)) {
                face::Pos<FaceType> p(&f, j);
                face::Pos<FaceType> startPos = p;
                assert(p.IsBorder());
                do {
                    assert(p.IsManifold());
                    p.F()->SetV();
                    //outline.push_back(Point2<ScalarType>(p.V()->P()));
                    outline.push_back(p.F()->WT(p.VInd()).P());
                    p.NextB();
                }
                while (p != startPos);
                outline2Vec.push_back(outline);
                outline.clear();
            }
        }
    }

    return outline2Vec.size();
}

/*
template <typename ScalarType = float>
void MeshOutlinesUV(Mesh& m, std::vector<std::vector<Point2<ScalarType>>> &outline2Vec)
{
    tri::UpdateFlags<Mesh>::FaceClearV(m);
    tri::UpdateFlags<Mesh>::VertexClearV(m);
    tri::UpdateTopology<Mesh>::FaceFace(m);

    outline2Vec.clear();
    std::vector<Point2<ScalarType>> outline;

    for (auto& f : m.face) {
        auto fptr = &f;
        for (int i = 0; i < 3; ++i) {
            if (!fptr->IsV() && face::IsBorder(*fptr, i)) {
                face::Pos<Mesh::FaceType> p(fptr, i);
                face::Pos<Mesh::FaceType> startPos = p;
                assert(p.IsBorder());
                do {
                    assert(p.IsManifold());
                    p.F()->SetV();
                    //outline.push_back(Point2<ScalarType>(p.V()->P()));
                    Point2d uv = p.F()->WT(p.VInd()).P();
                    outline.push_back(Point2<ScalarType>{ScalarType(uv[0]), ScalarType(uv[1])});
                    p.NextB();
                }
                while (p != startPos);
                outline2Vec.push_back(outline);
                outline.clear();
            }
        }
    }
}
*/

#endif // UV_H
