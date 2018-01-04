#ifndef UV_H
#define UV_H

#include <vcg/space/texcoord2.h>

#include <vector>
#include <memory>

#include <QImage>

//#include "mesh_graph.h"

struct TexCoordStorage {
    vcg::TexCoord2d tc[3];
};

// Save a copy of the original texture coordinates (this will be used to render the new texture)
template <class MeshType>
void StoreWedgeTexCoordAsAttribute(MeshType &m)
{
    //TODO assert MeshType has wedge tex coord
    typename MeshType::template PerFaceAttributeHandle<TexCoordStorage> WTCSh
            = tri::Allocator<MeshType>::template GetPerFaceAttribute<TexCoordStorage>(m, "WedgeTexCoordStorage");
    for (auto &f : m.face) {
        for (int i = 0; i < 3; ++i) {
            WTCSh[&f].tc[i].P() = f.WT(i).P();
            WTCSh[&f].tc[i].N() = f.WT(i).N();
        }
    }
}

// Computes per face connected component ids. Two face with the same id belong
// to the same connected component. Assumes the topology is already updated.
template <class MeshType>
std::size_t ComputePerFaceConnectedComponentIdAttribute(MeshType &m)
{
    typename MeshType::template PerFaceAttributeHandle<RegionID> CCIDh
            = tri::Allocator<MeshType>::template GetPerFaceAttribute<std::size_t>(m, "ConnectedComponentID");

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
static std::size_t ConvertTextureBoundaryToOutline2Vec(MeshType &m, std::vector<std::vector<Point2<ScalarType>>> &outline2Vec)
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

#endif // UV_H

