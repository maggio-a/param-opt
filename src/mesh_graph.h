#ifndef MESH_GRAPH_H
#define MESH_GRAPH_H

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>

#include <QImage>


#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>
#include <vcg/math/histogram.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/update/color.h>

#include "mesh.h"
#include "gl_util.h"
#include "distortion_pos.h"
#include "uv.h"

/// TODO cache border and area state
struct FaceGroup {
    Mesh& mesh;
    const RegionID id;
    std::vector<Mesh::FacePointer> fpVec;
    std::unordered_set<std::shared_ptr<FaceGroup>> adj;

    int numMerges;

    float minMappedFaceValue;
    float maxMappedFaceValue;

    float borderUV;
    bool borderChanged;

    FaceGroup(Mesh& m, const RegionID id_) : mesh{m}, id{id_}, fpVec{}, adj{}, numMerges{0},
        minMappedFaceValue{-1}, maxMappedFaceValue{-1}, borderUV{0}, borderChanged{false}
    {
    }

    void AddFace(const Mesh::FacePointer fptr)
    {
        fpVec.push_back(fptr);
        borderChanged = true;
    }

    float AreaUV() const
    {
        float areaUV = 0;
        for (auto fptr : fpVec) areaUV += std::abs(tri::Distortion<Mesh,true>::AreaUV(fptr));
        return areaUV;
    }

    float Area3D() const
    {
        float area3D = 0;
        for (auto fptr : fpVec) area3D += tri::Distortion<Mesh,true>::Area3D(fptr);
        return area3D;
    }

    float BorderUV()
    {
        if (borderChanged) UpdateBorder();
        return borderUV;
    }

    vcg::Box2f UVBox() const
    {
        vcg::Box2f box;
        for (auto fptr : fpVec) {
            box.Add(fptr->WT(0).P());
            box.Add(fptr->WT(1).P());
            box.Add(fptr->WT(2).P());
        }
        return box;
    }

    void NotifyParameterizationChange() { borderChanged = true; }

    Mesh::FacePointer Fp() { assert(!fpVec.empty()); return fpVec[0]; }

    std::size_t FN() const { return fpVec.size(); }
    std::size_t NumAdj() const { return adj.size(); }

    template <typename VP>
    void MapDistortionFromTexCoord(DistortionWedge::DistType distortionType, float areaScale, float edgeScale, const VP& vp)
    {
        minMappedFaceValue = std::numeric_limits<float>::max();
        maxMappedFaceValue = std::numeric_limits<float>::lowest();
        for (auto fptr : fpVec) {
            if (distortionType == DistortionWedge::AreaDist) {

                float areaUV = DistortionPos<Mesh>::AreaUV(fptr) * areaScale;
                float area3D = DistortionPos<Mesh>::Area3D(fptr, vp);
                assert(area3D > 0);
                float diff = (areaUV - area3D) / area3D;
                assert(!math::IsNAN(diff));
                fptr->Q() = diff;
                //fptr->Q() = DistortionWedge::AreaDistortion(fptr, areaScale);
            }
            else if (distortionType == DistortionWedge::AngleDist) {
                fptr->Q() = DistortionPos<Mesh>::AngleDistortion(fptr, vp);
            }
            else if (distortionType == DistortionWedge::EdgeDist) {
                fptr->Q() = (DistortionPos<Mesh>::EdgeDistortion(fptr, 0, edgeScale, vp) +
                             DistortionPos<Mesh>::EdgeDistortion(fptr, 1, edgeScale, vp) +
                             DistortionPos<Mesh>::EdgeDistortion(fptr, 2, edgeScale, vp)) / 3.0f;
            }
            else assert(0 && "FaceGroup::MapDistortion");
            minMappedFaceValue = std::min(minMappedFaceValue, fptr->Q());
            maxMappedFaceValue = std::max(maxMappedFaceValue, fptr->Q());
        }
    }

    void MapDistortion(DistortionWedge::DistType distortionType, float areaScale, float edgeScale)
    {
        minMappedFaceValue = std::numeric_limits<float>::max();
        maxMappedFaceValue = std::numeric_limits<float>::lowest();
        for (auto fptr : fpVec) {
            if (distortionType == DistortionWedge::AreaDist) {

                float areaUV = DistortionWedge::AreaUV(fptr) * areaScale;
                float area3D = DistortionWedge::Area3D(fptr);
                assert(area3D > 0);
                float diff = (areaUV - area3D) / area3D;
                assert(!math::IsNAN(diff));
                fptr->Q() = diff;
                //fptr->Q() = DistortionWedge::AreaDistortion(fptr, areaScale);
            }
            else if (distortionType == DistortionWedge::AngleDist) {
                fptr->Q() = DistortionWedge::AngleDistortion(fptr);
            }
            else if (distortionType == DistortionWedge::EdgeDist) {
                fptr->Q() = (DistortionWedge::EdgeDistortion(fptr, 0, edgeScale) +
                             DistortionWedge::EdgeDistortion(fptr, 1, edgeScale) +
                             DistortionWedge::EdgeDistortion(fptr, 2, edgeScale)) / 3.0f;
            }
            else assert(0 && "FaceGroup::MapDistortion");
            minMappedFaceValue = std::min(minMappedFaceValue, fptr->Q());
            maxMappedFaceValue = std::max(maxMappedFaceValue, fptr->Q());
        }
    }

    void UpdateBorder()
    {
        auto CCIDh = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(mesh, "ConnectedComponentID");
        assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(mesh, CCIDh));

        borderUV = 0.0f;
        for (auto fptr : fpVec) {
            for (int i = 0; i < 3; ++i) {
                if (face::IsBorder(*fptr, i) || CCIDh[fptr] != CCIDh[fptr->FFp(i)]) {
                    borderUV += DistortionWedge::EdgeLenghtUV(fptr, i);
                }
            }
        }
        borderChanged = false;
    }

};

template <typename MeshType, typename WedgeTextureMapper>
void CopyFaceGroupIntoMesh(MeshType &m, const struct FaceGroup& fg, std::unordered_map<Mesh::VertexPointer, typename MeshType::VertexPointer> &vpmap, const WedgeTextureMapper& wtmapper)
{
    m.Clear();
    vpmap.clear();
    vpmap.reserve(fg.FN() * 3);

    std::size_t vn = 0;
    for (auto fptr : fg.fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (vpmap.count(fptr->V(i)) == 0) {
                vn++;
                vpmap[fptr->V(i)] = nullptr;
            }
        }
    }
    auto mvi = tri::Allocator<MeshType>::AddVertices(m, vn);
    auto mfi = tri::Allocator<MeshType>::AddFaces(m, fg.FN());

    for (auto fptr : fg.fpVec) {
        typename MeshType::FacePointer mfp = &*mfi++;
        for (int i = 0; i < 3; ++i) {
            Mesh::VertexPointer vp = fptr->V(i);
            typename MeshType::VertexPointer& mvp = vpmap[vp];
            if (mvp == nullptr) {
                mvp = &*mvi++;
                mvp->P() = vp->P();
            }
            mfp->V(i) = mvp;
            mfp->WT(i) = wtmapper(fptr, i);
        }
    }
}


struct MeshGraph {

    Mesh& mesh;

    std::unordered_map<std::size_t, std::shared_ptr<FaceGroup>> charts;
    TextureObjectHandle textureObject;

    MeshGraph(Mesh& m) : mesh{m} {}

    void MapDistortionFromTexCoord(DistortionWedge::DistType distortionType)
    {
        auto WTCSh = tri::Allocator<Mesh>::FindPerFaceAttribute<TexCoordStorage>(mesh, "WedgeTexCoordStorage");
        assert(tri::Allocator<Mesh>::IsValidHandle<TexCoordStorage>(mesh, WTCSh));

        float areaScale;
        float edgeScale;
        auto VertexPosMapper = [&WTCSh] (Mesh::ConstFacePointer fptr, int i) {
            Point2f uv = WTCSh[fptr].tc[i].P(); return Point3f{uv[0], uv[1], 0};
        };

        DistortionPos<Mesh>::MeshScalingFactor(mesh, areaScale, edgeScale, VertexPosMapper);
        for (auto& c : charts) c.second->MapDistortionFromTexCoord(distortionType, areaScale, edgeScale, VertexPosMapper);
        Distribution<float> qd;
        tri::Stat<Mesh>::ComputePerFaceQualityDistribution(mesh, qd);
        std::pair<float, float> range = DistortionRange();
        tri::UpdateColor<Mesh>::PerFaceQualityRamp(mesh, range.first, range.second);
   }

    void MapDistortion(DistortionWedge::DistType distortionType)
    {
        float areaScale;
        float edgeScale;
        DistortionWedge::MeshScalingFactor(mesh, areaScale, edgeScale);
        for (auto& c : charts) c.second->MapDistortion(distortionType, areaScale, edgeScale);
        Distribution<float> qd;
        tri::Stat<Mesh>::ComputePerFaceQualityDistribution(mesh, qd);
        std::pair<float, float> range = DistortionRange();
        tri::UpdateColor<Mesh>::PerFaceQualityRamp(mesh, range.first, range.second);
   }

    std::pair<float,float> DistortionRange() const
    {
        std::pair<float,float> range = std::make_pair(std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest());
        for (const auto& c : charts) {
            range.first = std::min(c.second->minMappedFaceValue, range.first);
            range.second = std::max(c.second->maxMappedFaceValue, range.second);
        }
        return range;
    }

    std::shared_ptr<FaceGroup> GetChart(std::size_t i)
    {
        assert(charts.find(i) != charts.end() && "Chart does not exist");
        return charts[i];
    }

    std::shared_ptr<FaceGroup> GetChart_Insert(std::size_t i)
    {
        if (charts.find(i) == charts.end()) charts.insert(std::make_pair(i, std::make_shared<FaceGroup>(mesh, i)));
        return charts[i];
    }

    std::size_t Count() const
    {
        return charts.size();
    }

    int MergeCount() const
    {
        int n = 0;
        for (const auto& c : charts) n += c.second->numMerges;
        return n;

    }

    float Area3D() const
    {
        float area3D = 0.0f;
        for (const auto& c : charts) area3D += c.second->Area3D();
        return area3D;
    }

    float AreaUV() const
    {
        float areaUV = 0.0f;
        for (const auto& c : charts) areaUV += c.second->AreaUV();
        return areaUV;
    }

    /// TODO remove useless parameters
    float BorderUV(float *meshBorderLengthUV = nullptr, float *seamLengthUV = nullptr)
    {
        float borderUV = 0.0f;
        for (const auto& c : charts) borderUV += c.second->BorderUV();
        return borderUV;
    }
};

static void PrintParameterizationInfo(std::shared_ptr<MeshGraph> pdata)
{
    std::cout << pdata->charts.size() << " " << pdata->Area3D() << " "
              << pdata->AreaUV() << " " << pdata->BorderUV() << std::endl;
}

template <class MeshType>
std::shared_ptr<MeshGraph> ComputeParameterizationGraph(
        MeshType &m, TextureObjectHandle textureObject, float *uvMeshBorder = nullptr)
{
    std::size_t numRegions = ComputePerFaceConnectedComponentIdAttribute<MeshType>(m);

    std::shared_ptr<MeshGraph> paramData = std::make_shared<MeshGraph>(m);
    paramData->textureObject = textureObject;
    paramData->charts.reserve(numRegions);
    auto CCIDh = tri::Allocator<MeshType>::template GetPerFaceAttribute<std::size_t>(m, "ConnectedComponentID");

    auto ICCIDh = tri::Allocator<MeshType>::template GetPerFaceAttribute<std::size_t>(m, "InitialConnectedComponentID");
    ICCIDh._handle->data.assign(CCIDh._handle->data.begin(), CCIDh._handle->data.end());

    // build parameterization graph
    tri::UpdateTopology<Mesh>::FaceFace(m);
    for (auto &f : m.face) {
        std::size_t regionId = CCIDh[&f];
        paramData->GetChart_Insert(regionId)->AddFace(&f);
        // TODO this may be refactored into AddFace
        for (int i = 0; i < f.VN(); ++i) {
            std::size_t adjId = CCIDh[f.FFp(i)];
            if (regionId != adjId) {
                (paramData->GetChart_Insert(regionId)->adj).insert(paramData->GetChart_Insert(adjId));
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

    int Collapse_OK               = 0;
    int Collapse_ERR_DISCONNECTED = 1;
    int Collapse_ERR_UNFEASIBLE   = 2;

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


    /// TODO replace with a proper interface to the graph instance
    std::shared_ptr<MeshGraph> Graph()
    {
        return g;
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

    /// TODO this is not robust if a client tries to collapse an arbitrary edge rather than the one returned
    /// by *NextEdge() since feasibility is lazily evaluated (unfeasible edges may be present in the edge set)
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
    /// TODO reorder interface properly
private:
    void Merge(ChartHandle c1, ChartHandle c2)
    {
        assert(edges.find(Edge{c1,c2}) != edges.end());

        auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(g->mesh, "ConnectedComponentID");
        assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(g->mesh, CHIDh));

        // Update ID and add faces
        for (auto fp : c2->fpVec) {
            CHIDh[fp] = c1->id;
            c1->AddFace(fp);
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
                edges.erase(e1); // delete the edge if it already exists so it can be reweighted
                auto weighted = std::make_pair(e1, wfct(e1));
                edges.insert(weighted);
                queue.push(weighted);
            }
        }

        c1->numMerges += c2->numMerges + 1;

        g->charts.erase(c2->id);
    }
public:

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
                    // don't bother checking feasibility of each merge since we already know that the union of all the charts is
                    Merge(regions[entry.first], regions[id]);
                    mergeCount++;
                }
            } else {
                //tri::io::ExporterOBJ<PMesh>::Save(probe, "fail.obj", tri::io::Mask::IOM_WEDGTEXCOORD);
            }
        }

        return mergeCount;
    }

    /*
     * Test if a range of chart can be merged together
     * If the range can be merged together, it returns a queue defining a sequence of merges that result in the union of the chart range.
     * The merges can be performed by extracting the first chart, and iteratively merging it with the following charts in the queue.
     * If the range cannot be merged (either because it is disconnected or unfeasible) it returns a proper error code and an empty queue.
     * */
    template <typename ChartInputIterator>
    std::pair<int,std::queue<ChartHandle>> CollapseAllowed(ChartInputIterator first, ChartInputIterator last)
    {
        std::unordered_set<ChartHandle> chartSet, testSet;
        while (first != last) {
            if (*first) {
                chartSet.insert(*first);
                testSet.insert(*first);
            }
            first++;
        }

        // if chartset is empty or 1 chart it can always be merged
        if (chartSet.size() < 2) {
            auto p = std::make_pair(Collapse_OK, std::queue<ChartHandle>{});
            for (auto chart : chartSet) p.second.push(chart);
            return p;
        }

        // first validate the chart set by checking if an equivalent tree exists in the graph
        std::queue<ChartHandle> q;
        std::queue<ChartHandle> mergeQueue;
        q.push(*chartSet.begin());
        testSet.erase(q.front());
        while (q.size() > 0) {
            auto& chart = q.front();
            for (auto& c : chart->adj) {
                if (testSet.count(c) > 0) {
                    q.push(c);
                    testSet.erase(c);
                }
            }
            mergeQueue.push(chart);
            q.pop();
        }

        if (testSet.size() > 0) {
            return std::make_pair(Collapse_ERR_DISCONNECTED, std::queue<ChartHandle>{});
        }

        // build the test mesh, and merge the charts if the test mesh is feasible
        PMesh probe;
        std::vector<std::vector<Mesh::FacePointer>* > fpVecp;
        for (auto& chart : chartSet) {
            fpVecp.push_back(&(chart->fpVec));
        }
        BuildPMeshFromFacePointers<Mesh>(probe, fpVecp);
        if (Parameterizable(probe)) {
            return std::make_pair(Collapse_OK, mergeQueue);
        } else {
            //tri::io::ExporterOBJ<PMesh>::Save(probe, "fail.obj", tri::io::Mask::IOM_WEDGTEXCOORD);
            return std::make_pair(Collapse_ERR_UNFEASIBLE, std::queue<ChartHandle>{});
        }
    }

    // Collapse multiple charts in one go
    // if the merge cannot happen (either due to infeasibility of the resulting surface or because the charts are
    // not contiguous) it returns a proper error code
    template <typename ChartInputIterator>
    std::pair<int,ChartHandle> Collapse(ChartInputIterator first, ChartInputIterator last)
    {
        auto p = CollapseAllowed(first, last);
        if (p.first == Collapse_OK) {
            std::cout << "Merging charts" << std::endl;
            std::queue<ChartHandle>& mergeQueue = p.second;
            ChartHandle chart = mergeQueue.front();
            mergeQueue.pop();
            while (mergeQueue.size() > 0) {
                Merge(chart, mergeQueue.front());
                mergeQueue.pop();
            }
            return std::make_pair(Collapse_OK, chart);

        } else return std::make_pair(p.first, nullptr);
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
