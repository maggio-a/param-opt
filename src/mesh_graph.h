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
#include <vcg/complex/algorithms/hole.h>

#include "mesh.h"
#include "gl_utils.h"
#include "math_utils.h"
#include "metric.h"
#include "uv.h"


class MeshGraph;
class FaceGroup;

using ChartHandle = std::shared_ptr<FaceGroup>;

void PrintParameterizationInfo(std::shared_ptr<MeshGraph> pdata);

/*
 * Constructs the MeshGraph for the given mesh, returning a smart pointer as the handle to the object
 */
template <class MeshType>
std::shared_ptr<MeshGraph> ComputeParameterizationGraph(MeshType &m, TextureObjectHandle textureObject, float *uvMeshBorder = nullptr);

/*
 * Constructs a mesh from a FaceGroup. If sanitize == true closes the holes and splits non manifold vertices.
 * If wtcsattr == false does not copy the WedgeTexCoordStorage attribute into the new mesh. Extra faces added
 * in the hole filling process get an isometric parameterization scaled according to the chart attributes.
 * The added faces are visited.
 */
template <typename MeshType>
void CopyFaceGroupIntoMesh(MeshType &m, FaceGroup& fg, bool sanitize, bool wtcsattr);

/* Computes the texture outlines of a given chart */
template <typename ScalarType = float>
void ChartOutlinesUV(Mesh& m, ChartHandle chart, std::vector<std::vector<Point2<ScalarType>>> &outline2Vec);

/*
 * FaceGroup class
 *
 * Used to store a mesh chart as an array of Face pointers
 */
struct FaceGroup {
    Mesh& mesh;
    const RegionID id;
    std::vector<Mesh::FacePointer> fpVec;
    std::unordered_set<std::shared_ptr<FaceGroup>> adj;

    int numMerges;

    float minMappedFaceValue;
    float maxMappedFaceValue;

    mutable float borderUV;
    mutable bool borderChanged;

    FaceGroup(Mesh& m, const RegionID id_);

    void AddFace(const Mesh::FacePointer fptr);
    void ParameterizationChanged();

    Mesh::FacePointer Fp();

    std::size_t FN() const;
    std::size_t NumAdj() const;

    float AreaUV() const;
    float Area3D() const;
    float BorderUV() const;
    vcg::Box2d UVBox() const;


    void MapDistortion(DistortionMetric::Type type, ParameterizationGeometry geometry, double areaScale);
    void UpdateBorder() const;

};


/*
 * MeshGraph class
 *
 * The graph is actually stored as an associative array mapping each Region id to the relative FaceGroup, the adjacencies
 * are recorded inside each FaceGroup
 */
struct MeshGraph {

    Mesh& mesh;

    std::unordered_map<std::size_t, std::shared_ptr<FaceGroup>> charts;
    TextureObjectHandle textureObject;

    MeshGraph(Mesh& m);

    /* Map distortion values to the mesh graph, this essentially computes the distortion
     * for each face of the mesh, and then computes the minmax pair at each chart
     * param 'type' defines the type of distortion computed
     * param 'geometry' is one of
     *   ParameterizationGeometry::Model compute distortion wrt to the actual mesh geometry
     *   ParameterizationGeometry::Texture compute distortion wrt the original texture coordinates
     */
    void MapDistortion(DistortionMetric::Type type, ParameterizationGeometry geometry);

    /* compute the minmax distortion of the graph */
    std::pair<float,float> DistortionRange() const;

    /* Retrieve region i (assert if not found) */
    std::shared_ptr<FaceGroup> GetChart(RegionID i);

    /* Retrieve region i, creating if it is not found */
    std::shared_ptr<FaceGroup> GetChart_Insert(RegionID i);

    /* Number of regions */
    std::size_t Count() const;

    /* Number of merges performed on the graph (is it used?) */
    int MergeCount() const;

    float Area3D() const;
    float AreaUV() const;

    /// TODO remove useless parameters
    float BorderUV(float *meshBorderLengthUV = nullptr, float *seamLengthUV = nullptr);
};


class EdgeWeightFunction;

class GraphManager {

public:

    using ChartHandle = std::shared_ptr<FaceGroup>;

    static int Collapse_OK;
    static int Collapse_ERR_DISCONNECTED;
    static int Collapse_ERR_UNFEASIBLE;

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

    struct EdgeWeightComparator {
        int operator()(const std::pair<Edge,double>& e1, const std::pair<Edge,double>& e2) const
        {
            return e1.second > e2.second;
        }
    };

    GraphManager(std::shared_ptr<MeshGraph> gptr, std::unique_ptr<EdgeWeightFunction>&& wfct);

    /// TODO replace with a proper interface to the graph instance
    std::shared_ptr<MeshGraph> Graph();

    bool Valid(std::pair<Edge,double> weightedEdge);

    const std::pair<Edge,double>& PeekNextEdge();

    void RemoveNextEdge();

    bool HasNextEdge();

    /// TODO this is not robust if a client tries to collapse an arbitrary edge rather than the one returned
    /// by *NextEdge() since feasibility is lazily evaluated (unfeasible edges may be present in the edge set)
    ChartHandle Collapse(const Edge& e);

    // Splits the chart into its initial components, returns the handles to the created charts into splitCharts
    void Split(const RegionID id, std::vector<ChartHandle> &splitCharts);

    int CloseMacroRegions(std::size_t minRegionSize);

    /*
     * Test if a range of charts can be merged together
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

    /*
     * Collapse multiple charts in one go
     * if the merge cannot happen (either due to infeasibility of the resulting surface or because the charts are
     * not contiguous) it returns a proper error code
     */
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

    std::unique_ptr<EdgeWeightFunction> wfct;

    bool AddEdge(ChartHandle c1, ChartHandle c2, bool replace = false);
    void Merge(ChartHandle c1, ChartHandle c2);
};

struct EdgeWeightFunction {
    virtual ~EdgeWeightFunction() {}

    virtual double operator()(GraphManager::Edge& e) const = 0;
};

struct FaceSizeWeightedShared3DBorder : EdgeWeightFunction {

    static std::string Name() { return "FaceSizeWeighted3DBorder"; }

    Mesh& mesh;

    FaceSizeWeightedShared3DBorder(Mesh& m) : mesh{m}{}

    // The weight of an edge is the number of faces of the smallest chart divided by the fraction of the
    // total chart border shared with the other chart. This weight should prioritize smaller charts, and
    // among smaller charts the ones that, when merged, minimize the residual texture border
    double operator()(GraphManager::Edge& e) const
    {
        // First evaluate shared border
        auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<std::size_t>(mesh, "ConnectedComponentID");
        assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(mesh, CHIDh));

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

        double w = double(std::min(e.a->FN(), e.b->FN())) / (sharedBorder/totalBorder); // The smaller the shared fraction, the larger the weight
        assert(std::isfinite(w));
        return w;
    }
};




// Template functions implementation
// =================================

template <class MeshType>
std::shared_ptr<MeshGraph> ComputeParameterizationGraph(MeshType &m, TextureObjectHandle textureObject, float *uvMeshBorder)
{
    std::size_t numRegions = ComputePerFaceConnectedComponentIdAttribute<MeshType>(m);

    std::shared_ptr<MeshGraph> paramData = std::make_shared<MeshGraph>(m);
    paramData->textureObject = textureObject;
    paramData->charts.reserve(numRegions);
    auto CCIDh = tri::Allocator<MeshType>::template GetPerFaceAttribute<RegionID>(m, "ConnectedComponentID");

    auto ICCIDh = tri::Allocator<MeshType>::template GetPerFaceAttribute<RegionID>(m, "InitialConnectedComponentID");
    ICCIDh._handle->data.assign(CCIDh._handle->data.begin(), CCIDh._handle->data.end());

    // build parameterization graph
    tri::UpdateTopology<Mesh>::FaceFace(m);
    for (auto &f : m.face) {
        RegionID regionId = CCIDh[&f];
        paramData->GetChart_Insert(regionId)->AddFace(&f);
        // TODO this may be refactored into AddFace
        for (int i = 0; i < f.VN(); ++i) {
            RegionID adjId = CCIDh[f.FFp(i)];
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


template <typename MeshType>
void CopyFaceGroupIntoMesh(MeshType &m, FaceGroup& fg, bool sanitize, bool wtcsattr)
{
    m.Clear();
    std::unordered_map<Mesh::VertexPointer, typename MeshType::VertexPointer> vpmap;
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
            mfp->WT(i) = mvp->T();
        }
    }

    // face number before sanitizing mesh
    int startFN = m.FN();

    if (sanitize) {
        tri::UpdateTopology<MeshType>::FaceFace(m);

        int splitCount = tri::Clean<MeshType>::SplitNonManifoldVertex(m, 0);
        if (splitCount > 0) {
            std::cout << "Mesh was not vertex-manifold, " << splitCount << " vertices split" << std::endl;
        }

        std::vector<double> vTotalBorderLength;
        std::vector<std::vector<std::size_t>> vBorderFaces;

        tri::UpdateFlags<MeshType>::FaceClearV(m);
        for (auto& f : m.face) {
            for (int i = 0; i < 3; ++i) {
                if (!f.IsV() && face::IsBorder(f, i)) {
                    double totalBorderLength = 0;
                    std::vector<std::size_t> borderFaces;

                    face::Pos<typename MeshType::FaceType> p(&f, i);
                    face::Pos<typename MeshType::FaceType> startPos = p;
                    assert(p.IsBorder());
                    do {
                        assert(p.IsManifold());
                        p.F()->SetV();
                        borderFaces.push_back(tri::Index(m, p.F()));
                        totalBorderLength += EdgeLength(*p.F(), p.VInd());
                        p.NextB();
                    } while (p != startPos);
                    vTotalBorderLength.push_back(totalBorderLength);
                    vBorderFaces.push_back(borderFaces);
                }
            }
        }

        tri::UpdateFlags<MeshType>::FaceClearS(m);
        /// TODO in case multiple borders have the same size, it could be chosen randomly
        assert(vBorderFaces.size() > 0 && "Mesh has no boundaries");
        // select longest border (this is the only one that remains), while the others are closed
        if (vBorderFaces.size() > 1) {
            std::size_t k = std::distance(vTotalBorderLength.begin(), std::max_element(vTotalBorderLength.begin(), vTotalBorderLength.end()));

            for (std::size_t i = 0; i < vBorderFaces.size(); ++i) {
                if (i == k) continue;
                for (auto j : vBorderFaces[i]) m.face[j].SetS();
            }

            // close everything that is selected
            tri::Hole<MeshType>::template EarCuttingFill<tri::MinimumWeightEar<MeshType>>(m, m.FN(), true);

            /*
            // retrieve the hole size parameter, since all holes except the peripheral border, just use the border length minus 1
            for (std::size_t i = 0; i < vBorderVertices.size(); ++i) {
                if (i == k) continue;
                assert(vBorderVertices[k].size() > vBorderVertices[i].size());  // otherwise cannot close all holes below threshold
            }

            // close holes
            int maxHoleSize = vBorderVertices[k].size() - 1 + 1; // there is an off by one in the hole size computed by the vcglib
            tri::Hole<MeshType>::template EarCuttingFill<tri::MinimumWeightEar<MeshType>>(m, maxHoleSize, false);
            */
        }
    }

    /// FIXME XXX selecting faces to mark them as new...
    tri::UpdateFlags<MeshType>::FaceClearS(m);
    for (auto& f: m.face) {
        if (int(tri::Index(m, f)) >= startFN) f.SetS();
    }

    if (wtcsattr) {
        auto WTCSh = tri::Allocator<Mesh>::FindPerFaceAttribute<TexCoordStorage>(fg.mesh, "WedgeTexCoordStorage");
        assert(tri::Allocator<Mesh>::IsValidHandle<TexCoordStorage>(fg.mesh, WTCSh));
        auto WTCShNew = tri::Allocator<MeshType>::template GetPerFaceAttribute<TexCoordStorage>(m, "WedgeTexCoordStorage");

        /* if some holes were closed, the isometric parameterization of the faces is scaled to match the ratio
         * of the FaceGroup copied */
        double scale = std::sqrt(fg.AreaUV() / fg.Area3D());
        for (int i = 0; i < m.FN(); ++i) {
            if (i < startFN) {
                WTCShNew[m.face[i]] = WTCSh[fg.fpVec[i]];
            }
            else { // faces added to close holes, use isometric parameterization as attribute for those
                typename MeshType::FaceType& f = m.face[i];
                Point2d u10, u20;
                LocalIsometry(f.P(1) - f.P(0), f.P(2) - f.P(0), u10, u20);
                u10 *= scale; u20 *= scale;
                Point2d u0{0, 0};
                assert((u10 ^ u20) > 0);
                WTCShNew[&f].tc[0].P() = u0;
                WTCShNew[&f].tc[1].P() = u0 + u10;
                WTCShNew[&f].tc[2].P() = u0 + u20;
            }
        }
    }
}


/// FIXME if the chart is an aggregate that has not been reparameterized this breaks...
template <typename ScalarType>
void ChartOutlinesUV(Mesh& m, ChartHandle chart, std::vector<std::vector<Point2<ScalarType>>> &outline2Vec)
{
    assert(chart->numMerges == 0);

    struct FaceAdjacency {
        bool changed[3] = {false, false, false};
        Mesh::FacePointer FFp[3];
        char FFi[3];
    };

    auto CCIDh = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(m, CCIDh));

    std::unordered_map<Mesh::FacePointer,FaceAdjacency> savedAdj;
    for (auto fptr : chart->fpVec) {
        fptr->ClearV();
        for (int i = 0; i < 3; ++i) {
            if (CCIDh[fptr] != CCIDh[fptr->FFp(i)]) {
                savedAdj[fptr].changed[i] = true;
                savedAdj[fptr].FFp[i] = fptr->FFp(i);
                savedAdj[fptr].FFi[i] = fptr->FFi(i);
                fptr->FFp(i) = fptr;
                fptr->FFi(i) = i;
            }
        }
    }

    outline2Vec.clear();
    std::vector<Point2<ScalarType>> outline;

    for (auto fptr : chart->fpVec) {
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
                    outline.push_back(Point2<ScalarType>(uv[0], uv[1]));
                    p.NextB();
                }
                while (p != startPos);
                outline2Vec.push_back(outline);
                outline.clear();
            }
        }
    }

    for (auto& entry : savedAdj) {
        for (int i = 0; i < 3; ++i) {
            if (entry.second.changed[i]) {
                entry.first->FFp(i) = entry.second.FFp[i];
                entry.first->FFi(i) = entry.second.FFi[i];
            }
        }
    }

    std::size_t maxsz = 0;
    for (std::size_t i = 0; i < outline2Vec.size(); ++i) {
        maxsz = std::max(maxsz, outline2Vec[i].size());
    }
    if (maxsz == 0) {// fallback to uv box as outline
        vcg::Box2d box = chart->UVBox();
        std::vector<Point2<ScalarType>> outline;
        outline.push_back(Point2<ScalarType>(box.min.X(), box.min.Y()));
        outline.push_back(Point2<ScalarType>(box.max.X(), box.min.Y()));
        outline.push_back(Point2<ScalarType>(box.max.X(), box.max.Y()));
        outline.push_back(Point2<ScalarType>(box.min.X(), box.max.Y()));
        outline2Vec.push_back(outline);
    }

}


#endif // MESH_GRAPH_H
