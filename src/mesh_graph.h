#ifndef MESH_GRAPH_H
#define MESH_GRAPH_H

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>

#include <QImage>

#include "mesh.h"

struct FaceGroup {
    const RegionID id;
    std::vector<Mesh::FacePointer> fpVec;
    std::unordered_set<std::shared_ptr<FaceGroup>> adj;

    int numMerges;

    FaceGroup(const RegionID id_) : id{id_}, fpVec{}, adj{}, numMerges{0} {}

    void AddFace(const Mesh::FacePointer fptr, Mesh::PerFaceAttributeHandle<std::size_t>& CCIDh)
    {
        fpVec.push_back(fptr);
    }

    float AreaUV()
    {
        float areaUV = 0;
        for (auto fptr : fpVec) areaUV += std::abs(tri::Distortion<Mesh,true>::AreaUV(fptr));
        return areaUV;
    }

    float Area3D()
    {
        float area3D = 0;
        for (auto fptr : fpVec) area3D += tri::Distortion<Mesh,true>::Area3D(fptr);
        return area3D;
    }

    Mesh::FacePointer Fp() { assert(!fpVec.empty()); return fpVec[0]; }

    std::size_t FN() { return fpVec.size(); }
    std::size_t NumAdj() { return adj.size(); }
};

struct MeshGraph {

    Mesh& mesh;

    std::unordered_map<std::size_t, std::shared_ptr<FaceGroup>> charts;
    std::vector<std::shared_ptr<QImage>> textures;

    MeshGraph(Mesh& m) : mesh(m), charts{}, textures{} {}

    std::shared_ptr<FaceGroup> GetChart(std::size_t i)
    {
        if (charts.find(i) == charts.end()) charts.insert(std::make_pair(i, std::make_shared<FaceGroup>(i)));
        return charts[i];
    }

    std::size_t Count() {
        std::size_t sz = 0;
        for (auto c : charts) {
            sz += c.second->fpVec.size();
        }
        return sz;
    }

    float Area3D()
    {
        float area3D = 0.0f;
        for (auto c : charts) area3D += c.second->Area3D();
        return area3D;
    }

    float AreaUV()
    {
        float areaUV = 0.0f;
        for (auto c : charts) areaUV += c.second->AreaUV();
        return areaUV;
    }

    float BorderUV(float *meshBorderLengthUV = nullptr, float *seamLengthUV = nullptr)
    { std::cout << "WARNING: ParameterizationData::BorderUV() not implemented" << std::endl; return 0; }
};

/*
void ReparameterizeZeroAreaRegions(Mesh &m, ParameterizationData& pdata)
{
    float meshArea3D = pdata.Area3D();
    for (auto& entry : pdata.charts) {
        auto chart = entry.second;

        if (chart->AreaUV() > 0) continue;

        PMesh pm;
        for (auto fptr : chart->fpVec) {
            tri::Allocator<PMesh>::AddFace(pm, fptr->P(0), fptr->P(1), fptr->P(2));
        }

        tri::Clean<PMesh>::RemoveDuplicateVertex(pm);
        tri::Allocator<PMesh>::CompactEveryVector(pm);
        tri::UpdateTopology<PMesh>::FaceFace(pm);
        tri::UpdateSelection<PMesh>::Clear(pm);

        tri::PoissonSolver<PMesh> solver(pm);
        tri::UpdateBounding<PMesh>::Box(pm);
        if (!solver.IsFeasible()) {
            tri::io::ExporterOBJ<PMesh>::Save(pm, "debug-reparam.obj", tri::io::Mask::IOM_WEDGTEXCOORD);
            assert(0 && "Poisson solver unfeasible");
        }
        solver.Init();
        solver.FixDefaultVertices();
        solver.SolvePoisson();
        tri::UpdateTexture<PMesh>::WedgeTexFromVertexTex(pm);

        float uvArea = 0.0f;
        for (auto &pf : pm.face) {
            uvArea += tri::Distortion<PMesh,true>::AreaUV(&pf);
        }

        // attempt to make the uv area of the region somewhat proportional in uv space
        // to the surface area in 3D space
        float scale = std::sqrt(chart->Area3D() / (uvArea * meshArea3D));

        assert(scale > 0);

        for (std::size_t i = 0; i < pm.face.size(); ++i) {
            auto &pf = pm.face[i];
            auto &f = *(chart->fpVec[i]);
            for (int k = 0; k < f.VN(); ++k) {
                pf.WT(k).P() = pf.WT(k).P() * scale;
                f.WT(k) = pf.WT(k);
            }
        }

    }
} */


#endif // MESH_GRAPH_H
