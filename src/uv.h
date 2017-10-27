#ifndef UV_H
#define UV_H

#include <vcg/space/texcoord2.h>

#include <vector>
#include <memory>

#include <QImage>

#include "mesh_graph.h"

struct WedgeTexCoordStorage {
    vcg::TexCoord2f tc[3];
};

// Save a copy of the original texture coordinates (this will be used to render the new texture)
template <class MeshType>
void StoreWedgeTexCoordAsAttribute(MeshType &m)
{
    //TODO assert MeshType has wedge tex coord
    typename MeshType::template PerFaceAttributeHandle<WedgeTexCoordStorage> WTCSh
            = tri::Allocator<MeshType>::template GetPerFaceAttribute<WedgeTexCoordStorage>(m, "WedgeTexCoordStorage");
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

template <class MeshType>
std::shared_ptr<MeshGraph> ComputeParameterizationGraph(
        MeshType &m, std::vector<std::shared_ptr<QImage>> imgVec, float *uvMeshBorder = nullptr)
{
    std::size_t numRegions = ComputePerFaceConnectedComponentIdAttribute<MeshType>(m);

    std::shared_ptr<MeshGraph> paramData = std::make_shared<MeshGraph>(m);
    paramData->textures = imgVec;
    paramData->charts.reserve(numRegions);
    auto CCIDh = tri::Allocator<MeshType>::template GetPerFaceAttribute<std::size_t>(m, "ConnectedComponentID");

    auto ICCIDh = tri::Allocator<MeshType>::template GetPerFaceAttribute<std::size_t>(m, "InitialConnectedComponentID");
    ICCIDh._handle->data.assign(CCIDh._handle->data.begin(), CCIDh._handle->data.end());

    // build parameterization graph
    tri::UpdateTopology<Mesh>::FaceFace(m);
    for (auto &f : m.face) {
        std::size_t regionId = CCIDh[&f];
        paramData->GetChart(regionId)->AddFace(&f, CCIDh);
        // TODO this may be refactored into AddFace
        for (int i = 0; i < f.VN(); ++i) {
            std::size_t adjId = CCIDh[f.FFp(i)];
            if (regionId != adjId) {
                (paramData->GetChart(regionId)->adj).insert(paramData->GetChart(adjId));
            }
        }
    }

    // compute uv mesh border if required
    if (uvMeshBorder) {
        *uvMeshBorder = 0.0f;
        for (auto &f : m.face) {
            for (int i = 0; i < f.VN(); ++i) {
                if (face::IsBorder(f, i)) {
                   *uvMeshBorder += (f.cWT((i+1)%f.VN()).P() - f.cWT(i).P()).Norm();
                }
            }
        }
    }

    return paramData;
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

void PrintParameterizationInfo(std::shared_ptr<MeshGraph> pdata)
{
    std::cout << pdata->charts.size() << " " << pdata->Area3D() << " "
              << pdata->AreaUV() << " " << pdata->BorderUV() << std::endl;
}

#endif // UV_H

