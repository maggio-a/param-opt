#include "mesh.h"
#include "uv.h"
#include "metric.h"
#include "math_utils.h"
#include "mesh_attribute.h"
#include "mesh_graph.h"

using namespace vcg;

void StoreWedgeTexCoordAsAttribute(Mesh &m, GraphHandle graph)
{
    TextureObjectHandle textureObject = graph->textureObject;
    double uvRatio = textureObject->TextureWidth(0) / (double) textureObject->TextureHeight(0);

    auto WTCSh = GetWedgeTexCoordStorageAttribute(m);
    for (auto &f : m.face) {
        for (int i = 0; i < 3; ++i) {
            WTCSh[&f].tc[i].P() = f.WT(i).P();
            WTCSh[&f].tc[i].N() = f.WT(i).N();

            // Handle squashed texcoords. Note that this is the opposite of
            // what is done to parameterize zero uv area regions.
            if (uvRatio > 1)
                WTCSh[&f].tc[i].P().X() *= uvRatio;
            else if (uvRatio < 1)
                WTCSh[&f].tc[i].P().Y() /= uvRatio;
        }
    }
}

Box2d UVBox(const Mesh& m)
{
    Box2d uvbox;
    for(auto const& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            uvbox.Add(f.cWT(i).P());
        }
    }
    return uvbox;
}

Box2d UVBoxVertex(const Mesh& m)
{
    Box2d uvbox;
    for(auto const& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            uvbox.Add(f.cV(i)->T().P());
        }
    }
    return uvbox;
}


double ComputeParameterizationScaleInfo(Mesh& m)
{
    auto info = GetParameterizationScaleInfoAttribute(m);
    for (auto& f : m.face) {
        double surfaceAreaFace = DistortionMetric::Area3D(f);
        double paramAreaFace = std::abs(DistortionMetric::AreaUV(f));
        if (surfaceAreaFace > 0.0 && paramAreaFace > 0.0) {
            info().surfaceArea += surfaceAreaFace;
            info().parameterArea += paramAreaFace;
            info().numNonZero++;
        }
    }
    info().scale = std::sqrt(info().parameterArea / info().surfaceArea);
    return info().scale;
}

std::size_t ComputePerFaceConnectedComponentIdAttribute(Mesh &m)
{
    auto CCIDh = GetConnectedComponentIDAttribute(m);

    tri::UpdateFlags<Mesh>::FaceClearV(m);
    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(m);

    std::stack<Mesh::FacePointer> s;
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

void MarkSeamsAsFaux(Mesh& m)
{
    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(m);
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            if (face::IsBorder(f, i)) {
                f.SetF(i);
            }
        }
    }
}
