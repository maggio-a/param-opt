#ifndef MESH_GRAPH_H
#define MESH_GRAPH_H

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <algorithm>

#include <QImage>

#include "mesh.h"
#include "gl_utils.h"
#include "math_utils.h"
#include "metric.h"
#include "packing_utils.h"


struct MeshGraph;
struct FaceGroup;
struct EdgeWeightFunction;

using ChartHandle = std::shared_ptr<FaceGroup>;
using GraphHandle = std::shared_ptr<MeshGraph>;

void PrintParameterizationInfo(std::shared_ptr<MeshGraph> pdata);

/* Constructs the MeshGraph for the given mesh, returning a smart pointer as the
 * handle to the object */
std::shared_ptr<MeshGraph> ComputeParameterizationGraph(Mesh &m, TextureObjectHandle textureObject);

/* Constructs a mesh from a FaceGroup, the created mesh has the FaceIndex
 * attribute defined (see mesh_attribute.h) */
void CopyToMesh(FaceGroup& fg, Mesh& m);

void ExtractOutlines(GraphHandle graph, std::vector<Outline2f>& texOutlines, std::unordered_map<RegionID,std::size_t>& outlineIndex);

/* Computes the UV outline(s) of the given chart. If the chart has no outlines,
 * which can happen for some inputs on small closed components that are ignored
 * by the reparameterization procedure, it returns as outline the bounding box
 * of the chart texture coordinates.
 * NOTE: It assumes the face-face topology is computed according to the wedge
 * texture coordinates of the chart/mesh */
/// FIXME if the chart is an aggregate that has not been reparameterized this breaks...
void ChartOutlinesUV(Mesh& m, FaceGroup& chart, std::vector<std::vector<Point2d>> &outline2Vec);
void ChartOutlinesUV(Mesh& m, FaceGroup& chart, std::vector<std::vector<Point2f>> &outline2Vec);

/* Tries to remove topological noise in order to simplify the aggregation step.
 * This function potentially changes the mesh geometry, leaving to the calling
 * code the responsibility to recompute the mesh graph.
 * Returns true if the mesh is changed (therefore invalidating the graph object,
 * false otherwise. */
bool CleanSmallComponents(Mesh& m, GraphHandle graph, TextureObjectHandle texObj, double areaThreshold);

/* FaceGroup class
 * Used to store a mesh chart as an array of Face pointers */
struct FaceGroup {

    struct Hasher {
        std::size_t operator()(const ChartHandle& ch) const
        {
            return std::hash<RegionID>()(ch->id);
        }
    };

    struct Cache {
        double areaUV;
        double area3D;
        double borderUV;
        double border3D;
        vcg::Point3d weightedSumNormal;
    };

    void UpdateCache() const;

    Mesh& mesh;
    RegionID id;
    std::vector<Mesh::FacePointer> fpVec;
    std::unordered_set<ChartHandle, Hasher> adj;

    int numMerges;

    float minMappedFaceValue;
    float maxMappedFaceValue;

    double error;

    mutable bool dirty;
    mutable Cache cache;

    FaceGroup(Mesh& m, const RegionID id_);

    void Clear();
    void AddFace(const Mesh::FacePointer fptr);
    void ParameterizationChanged();
    Mesh::FacePointer Fp();

    vcg::Point3d AverageNormal() const;

    std::size_t FN() const;
    std::size_t NumAdj() const;
    double OriginalAreaUV() const;
    double AreaUV() const;
    double Area3D() const;
    double BorderUV() const;
    double Border3D() const;
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

    std::unordered_map<RegionID, std::shared_ptr<FaceGroup>> charts;
    TextureObjectHandle textureObject;

    MeshGraph(Mesh& m);
    ~MeshGraph();

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

    /* Retrieve region i (ensure_condition if not found) */
    std::shared_ptr<FaceGroup> GetChart(RegionID i);

    /* Retrieve region i, creating if it is not found */
    std::shared_ptr<FaceGroup> GetChart_Insert(RegionID i);

    /* Number of regions */
    std::size_t Count() const;

    /* Number of merges performed on the graph (is it used?) */
    int MergeCount() const;

    double Area3D() const;
    double AreaUV() const;

    double MappedFraction() const;

    double BorderUV() const;

    void ToEdgeMesh(MyMesh& em)
    {
        em.Clear();
        auto CCIDh = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(mesh, "ConnectedComponentID");
        ensure_condition(tri::Allocator<Mesh>::IsValidHandle<RegionID>(mesh, CCIDh));

        for (auto& f : mesh.face) {
            for (int i = 0; i < 3; ++i) {
                if (face::IsBorder(f, i) || (CCIDh[f] != CCIDh[f.FFp(i)])) {
                    tri::Allocator<MyMesh>::AddEdge(em, f.P0(i), f.P1(i));
                }
            }
        }
        tri::Clean<MyMesh>::RemoveDuplicateEdge(em);
        tri::Clean<MyMesh>::RemoveDuplicateVertex(em);
        tri::Allocator<MyMesh>::CompactEveryVector(em);
    }
};

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
            ensure_condition(a != b);
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

    int CloseMacroRegions(double areaThreshold);

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
        int fn = 0;
        for (auto& chart : chartSet)
            fn += chart->FN();
        std::vector<Mesh::FacePointer> fpv;
        fpv.reserve(fn);
        for (auto& chart : chartSet)
            fpv.insert(fpv.end(), chart->fpVec.begin(), chart->fpVec.end());

        Mesh probe;
        MeshFromFacePointers(fpv, probe);

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
            std::queue<ChartHandle>& mergeQueue = p.second;
            ChartHandle chart = mergeQueue.front();
            mergeQueue.pop();
            while (mergeQueue.size() > 0) {
                Merge(chart, mergeQueue.front());
                mergeQueue.pop();
            }
            return std::make_pair(Collapse_OK, chart);

        } else {
            return std::make_pair(p.first, nullptr);
        }
    }

    bool ExistsEdge(GraphManager::Edge e);

    double EdgeWeight(GraphManager::Edge e);

    void Debug_PrintGraphState();

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
    void RemoveEdge(ChartHandle c1, ChartHandle c2);
    void Merge(ChartHandle c1, ChartHandle c2);
};

struct EdgeWeightFunction {
    virtual ~EdgeWeightFunction() {}

    virtual double operator()(GraphManager::Edge& e) const = 0;
    virtual std::string Name() = 0;
};

struct W3D : EdgeWeightFunction {

    std::string Name() { return "w_3D"; }

    Mesh& mesh;

    W3D(Mesh& m) : mesh{m}{}

    double operator()(GraphManager::Edge& e) const
    {
        constexpr int FACE_THRESHOLD = 100000;

        bool islandEdge = ((e.a->NumAdj() == 1) || (e.b->NumAdj() == 1));

        if (islandEdge || ((e.a->FN() + e.b->FN()) < FACE_THRESHOLD)) {
            // First evaluate shared border
            auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<std::size_t>(mesh, "ConnectedComponentID");
            ensure_condition(tri::Allocator<Mesh>::IsValidHandle<RegionID>(mesh, CHIDh));

            double area_a = e.a->Area3D();
            double area_b = e.b->Area3D();
            ChartHandle chart1 = area_a < area_b ? e.a : e.b;
            ChartHandle chart2 = chart1 == e.a ? e.b : e.a;

            double totalBorder = 0.0f;
            double sharedBorder = 0.0f;
            for (auto fptr : chart1->fpVec) {
                for (int i = 0; i < fptr->VN(); ++i) {
                    auto ffpi = fptr->FFp(i);
                    auto adjId = CHIDh[ffpi];
                    if (adjId != CHIDh[fptr]) {
                        double edgeLength = EdgeLength(*fptr, i);
                        if (adjId == chart2->id)
                            sharedBorder += edgeLength;
                        totalBorder += edgeLength;
                    }
                }
            }

            double w = double(std::min(area_a, area_b)) * std::pow((sharedBorder/totalBorder), -1.0); // The smaller the shared fraction, the larger the weight
            //double w = double(std::min(area_a, area_b)) * (1.0 - (sharedBorder/totalBorder)); // The smaller the shared fraction, the larger the weight
            ensure_condition(std::isfinite(w));
            return w;
        } else {
            return Infinity();
        }

    }
};

struct W_Geometry3D : EdgeWeightFunction {

    std::string Name() { return "w_geometry_3D"; }

    Mesh& mesh;

    W_Geometry3D(Mesh& m) : mesh{m}{}

    double operator()(GraphManager::Edge& e) const
    {
        auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<std::size_t>(mesh, "ConnectedComponentID");
        ensure_condition(tri::Allocator<Mesh>::IsValidHandle<RegionID>(mesh, CHIDh));

        double border3D_contribution_a = e.a->Border3D();
        double border3D_contribution_b = e.b->Border3D();
        RegionID id_a = e.a->id;
        RegionID id_b = e.b->id;

        for (auto fptr : e.a->fpVec) {
            for (int i = 0; i < fptr->VN(); ++i) {
                auto ffpi = fptr->FFp(i);
                if (id_b == CHIDh[ffpi])
                    border3D_contribution_a -= EdgeLength(*fptr, i);
            }
        }

        for (auto fptr : e.b->fpVec) {
            for (int i = 0; i < fptr->VN(); ++i) {
                auto ffpi = fptr->FFp(i);
                if (id_a == CHIDh[ffpi])
                    border3D_contribution_b -= EdgeLength(*fptr, i);
            }
        }

        //ensure_condition(border3D >= 0);

        W3D w3d{mesh};
        double w3dw = w3d(e);

        vcg::Point3d n_a = e.a->AverageNormal();
        vcg::Point3d n_b = e.b->AverageNormal();

        //double w = (1 - (n_a.dot(n_b))) * border3D;
        //double w = (1 - (n_a.dot(n_b))) * std::min(border3D_contribution_a, border3D_contribution_b);
        double w = (1 - (n_a.dot(n_b))) * w3dw;

        return w;
    }
};



struct WFN : EdgeWeightFunction {

    //std::string Name() { return "FaceSizeWeighted3DBorder"; }
    std::string Name() { return "w_FN"; }

    Mesh& mesh;

    WFN(Mesh& m) : mesh{m}{}

    // The weight of an edge is the number of faces of the smallest chart divided by the fraction of the
    // total chart border shared with the other chart. This weight should prioritize smaller charts, and
    // among smaller charts the ones that, when merged, minimize the residual texture border
    double operator()(GraphManager::Edge& e) const
    {
        // First evaluate shared border
        auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<std::size_t>(mesh, "ConnectedComponentID");
        ensure_condition(tri::Allocator<Mesh>::IsValidHandle<RegionID>(mesh, CHIDh));

        ChartHandle chart1 = e.a->FN() < e.b->FN() ? e.a : e.b;
        ChartHandle chart2 = chart1 == e.a ? e.b : e.a;

        float totalBorder = 0.0f;
        float sharedBorder = 0.0f;
        for (auto fptr : chart1->fpVec) {
            for (int i = 0; i < fptr->VN(); ++i) {
                auto ffpi = fptr->FFp(i);
                auto adjId = CHIDh[ffpi];
                if (adjId != CHIDh[fptr]) {
                    float edgeLength = EdgeLength(*fptr, i);
                    if (adjId == chart2->id)
                        sharedBorder += edgeLength;
                    totalBorder += edgeLength;
                }
            }
        }

        double w = double(std::min(e.a->FN(), e.b->FN())) / (sharedBorder/totalBorder); // The smaller the shared fraction, the larger the weight
        ensure_condition(std::isfinite(w));
        return w;
    }
};

// if edge (a,b) with area uv a < b, then weight is (1-shared_uv_fraction(a)) * area_uv_a * abs(shared_fraction_a - shared_fraction_b)
struct WUV : EdgeWeightFunction {

    std::string Name() { return "w_UV"; }

    Mesh& mesh;

    WUV(Mesh& m) : mesh{m}{}

    double operator()(GraphManager::Edge& e) const
    {
        auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<std::size_t>(mesh, "ConnectedComponentID");
        ensure_condition(tri::Allocator<Mesh>::IsValidHandle<RegionID>(mesh, CHIDh));

        ChartHandle chart1 = e.a->FN() < e.b->FN() ? e.a : e.b;
        ChartHandle chart2 = chart1 == e.a ? e.b : e.a;

        double borderFrom1 = 0.0;
        double borderFrom2 = 0.0;
        for (auto fptr : chart1->fpVec) {
            for (int i = 0; i < fptr->VN(); ++i) {
                auto ffpi = fptr->FFp(i);
                auto adjId = CHIDh[ffpi];
                if (adjId == chart2->id) {
                    borderFrom1 += EdgeLengthUV(*fptr, i);
                    borderFrom2 += EdgeLengthUV(*ffpi, fptr->FFi(i));
                }
            }
        }

        double sharedFraction1 = borderFrom1 / chart1->BorderUV();
        double sharedFraction2 = borderFrom2 / chart2->BorderUV();

        double w1 = (1.0 - sharedFraction1); // if the border from 1 is fully shared with 2 the weight collapses to 0
        double w2 = chart1->AreaUV(); // area weighted, favor small regions in uv space
        double w3 = std::pow(borderFrom1 - borderFrom2, 2.0); // penalize merging regions with different resolution at the boundary in uv space

        double w = w1*w2*w3;
        ensure_condition(std::isfinite(w));
        return w;
    }
};

#endif // MESH_GRAPH_H
