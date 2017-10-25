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

class GraphManager {

public:

    struct Edge {
        std::shared_ptr<FaceGroup> a;
        std::shared_ptr<FaceGroup> b;

        Edge(std::shared_ptr<FaceGroup> e1, std::shared_ptr<FaceGroup> e2) : a{e1}, b{e2}
        {
            assert(a != b);
            if (b->id > a->id) std::swap(a, b);
        }

        bool operator==(const Edge& other) const noexcept
        {
            return a == other.a && b == other.b;
        }
    };

    struct EdgeHasher {
        std::size_t operator()(const Edge& e) const noexcept
        {
            std::size_t seed = 0;
            seed ^= std::hash<std::shared_ptr<FaceGroup>>()(e.a) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<std::shared_ptr<FaceGroup>>()(e.b) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };

    // TODO This will eventually become a template component
    struct EdgeWeightFunction {
        std::shared_ptr<MeshGraph> g;

        // The weight of an edge is the number of faces of the smallest chart divided by the fraction of the
        // total chart border shared with the other chart. This weight should prioritize smaller charts, and
        // among smaller charts the ones that, when merged, minimize the residual texture border
        double operator()(Edge& e) const
        {
            // First evaluate shared border
            auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<std::size_t>(g->mesh, "ConnectedComponentID");
            assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(g->mesh, CHIDh));

            ChartHandle chart1 = e.a->FN() < e.b->FN() ? e.a : e.b;
            ChartHandle chart2 = chart1 == e.a ? e.b : e.a;

            float totalBorder = 0.0f;
            float sharedBorder = 0.0f;
            for (auto fptr : chart1->fpVec) {
                for (int i = 0; i < fptr->VN(); ++i) {
                    auto ffpi = fptr->FFp(i);
                    auto adjId = CHIDh[ffpi];
                    if (adjId != CHIDh[fptr]) {
                        float edgeLength = DistortionWedge::EdgeLenght3D(fptr, i);
                        if (adjId == chart2->id)
                            sharedBorder += edgeLength;
                        totalBorder += edgeLength;
                    }
                }
            }

            return double(std::min(e.a->FN(), e.b->FN())) / (sharedBorder/totalBorder); // The smaller the shared fraction, the larger the weight
        }
    };

    using ChartHandle = std::shared_ptr<FaceGroup>;

    struct EdgeWeightComparator {
        int operator()(const std::pair<Edge,double>& e1, const std::pair<Edge,double>& e2) const
        {
            return e1.second > e2.second;
        }
    };

private:

    std::shared_ptr<MeshGraph> g;

    // Edge set (map), this is always consistent with the state of the mesh graph (that is, an edge exists in this
    // collection if and only if it exists in the mesh graph, and it is always mapped to the correct weight)
    std::unordered_map<Edge, double, EdgeHasher> edges;

    // Priority queue for the weighted edges. Elements in this queue may be no longer valid due to merge operations
    // (their weight might have changed, or the edge itself may no longer exist if it references a region that was
    // absorbed by another), so an edge in this queue is valid only if it exists in the collection _edges_ and the
    // weights are consistent
    std::priority_queue<std::pair<Edge,double>, std::vector<std::pair<Edge,double>>, EdgeWeightComparator> queue;

    EdgeWeightFunction wfct;

public:

    GraphManager(std::shared_ptr<MeshGraph> gptr) : g{gptr}, wfct{gptr}
    {
        tri::UpdateTopology<Mesh>::FaceFace(g->mesh);
        for (auto elem : g->charts) {
            ChartHandle c1 = elem.second;
            for (ChartHandle c2 : c1->adj) {
                Edge e{c1, c2};
                if (edges.find(e) == edges.end()) {
                    auto weightedEdge = std::make_pair(e, wfct(e));
                    edges.insert(weightedEdge);
                    queue.push(weightedEdge);
                }
            }
        }
    }

    bool Valid(std::pair<Edge,double> weightedEdge)
    {
        auto e = edges.find(weightedEdge.first);
        return e != edges.end() && e->second == weightedEdge.second;
    }

    const std::pair<Edge,double>& PeekNextEdge()
    {
        assert(HasNextEdge());
        return queue.top();
    }

    void RemoveNextEdge()
    {
        assert(HasNextEdge());
        queue.pop();
    }

    bool HasNextEdge()
    {
        constexpr double infinity = std::numeric_limits<double>::infinity();
        while (queue.size() > 0)
        {
            auto we = queue.top();
            if (Valid(we)) {
                if (we.second == infinity) { // The remaining edges cannot be collapsed due to infeasibility
                    return false;
                }
                else {
                    PMesh test;
                    Edge e = we.first;
                    BuildPMeshFromFacePointers<Mesh>(test, {&(e.a->fpVec), &(e.b->fpVec)});
                    if (Parameterizable(test)) {
                        return true;
                    }
                    else { // Cannot collapse the edge, reinsert with infinite weight
                        queue.pop();
                        we.second = infinity;
                        edges[we.first] = infinity;
                        queue.push(we);
                    }
                }
            }
            else { // Invalid edge, throw it away
                queue.pop();
            }
        }
        return false;
    }

    ChartHandle Collapse(const Edge& e)
    {
        ChartHandle c1, c2;
        if (e.a->FN() > e.b->FN()) {
            c1 = e.a; c2 = e.b;
        }
        else {
            c1 = e.b; c2 = e.a;
        }

        Merge(c1, c2);
        return c1;
    }

    // Merges c2 into c1
    void Merge(ChartHandle c1, ChartHandle c2)
    {
        assert(edges.find(Edge{c1,c2}) != edges.end());

        auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(g->mesh, "ConnectedComponentID");
        assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(g->mesh, CHIDh));

        // Update ID attribute
        for (auto fp : c2->fpVec) {
            CHIDh[fp] = c1->id;
            c1->AddFace(fp, CHIDh);
        }

        // Update adjacencies - Note that obsolete edges remain in the priority queue
        c1->adj.erase(c2);
        for (auto cn : c2->adj) {
            Edge e2{c2, cn};
            edges.erase(e2);
            if (cn != c1) { // c1 is now (if it wasn't already) adjacent to cn
                cn->adj.erase(c2);
                cn->adj.insert(c1);
                c1->adj.insert(cn);

                // Manage queue and edge map state
                Edge e1{c1, cn};
                edges.erase(e1); // delete the edge if it already exists as it
                auto weighted = std::make_pair(e1, wfct(e1));
                edges.insert(weighted);
                queue.push(weighted);
            }
        }

        c1->numMerges++;

        g->charts.erase(c2->id);
    }

    int CloseMacroRegions(std::size_t minRegionSize)
    {
        int mergeCount = 0;

        auto& regions = g->charts;

        std::unordered_map<RegionID, std::vector<RegionID>> mergeLists;
        std::unordered_map<RegionID, RegionID> invertedIndex;

        for (const auto& entry : regions) {
            auto chart = entry.second;
            if (invertedIndex.count(chart->id) == 1) continue; // skip if this region is already going to be merged to something else
            for (auto& adjRegion : chart->adj)
                if (adjRegion->NumAdj() == 1 && adjRegion->FN() < minRegionSize) {
                    assert(invertedIndex.count(adjRegion->id) == 0);
                    mergeLists[chart->id].push_back(adjRegion->id);
                    invertedIndex[adjRegion->id] = chart->id;
                }
        }

        for (auto& entry : mergeLists) {
            PMesh probe;
            std::vector<std::vector<Mesh::FacePointer>* > fpVecp;
            fpVecp.reserve(entry.second.size()+1);
            fpVecp.push_back(&(regions[entry.first]->fpVec));
            for (auto id : entry.second) {
                fpVecp.push_back(&(regions[id]->fpVec));
            }
            BuildPMeshFromFacePointers<Mesh>(probe, fpVecp);
            if (Parameterizable(probe)) {
                std::cout << "Merging " << entry.second.size() << " islands from macro region " << entry.first << std::endl;
                for (auto id : entry.second) {
                    Merge(regions[entry.first], regions[id]);
                    mergeCount++;
                }
            } else {
                //tri::io::ExporterOBJ<PMesh>::Save(probe, "fail.obj", tri::io::Mask::IOM_WEDGTEXCOORD);
            }
        }

        return mergeCount;
    }
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
