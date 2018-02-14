#include "mesh_graph.h"

#include "mesh.h"
#include "gl_utils.h"
#include "math_utils.h"
#include "metric.h"
#include "uv.h"

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>
#include <vcg/math/histogram.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/update/color.h>

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>

#include <QImage>

int GraphManager::Collapse_OK = 0;
int GraphManager::Collapse_ERR_DISCONNECTED = 1;
int GraphManager::Collapse_ERR_UNFEASIBLE = 2;

void PrintParameterizationInfo(std::shared_ptr<MeshGraph> pdata)
{
    std::cout << "Parameterization Info (Charts, Area 3D, Area UV, Border UV) " << std::endl;
    std::cout << pdata->charts.size() << " " << pdata->Area3D() << " "
              << pdata->AreaUV() << " " << pdata->BorderUV() << std::endl;
}


// FaceGroup class implementation
// ==============================

FaceGroup::FaceGroup(Mesh& m, const RegionID id_)
    : mesh{m},
      id{id_},
      fpVec{},
      adj{},
      numMerges{0},
      minMappedFaceValue{-1},
      maxMappedFaceValue{-1},
      borderUV{0},
      borderChanged{false}
{
}

void FaceGroup::AddFace(const Mesh::FacePointer fptr)
{
    fpVec.push_back(fptr);
    ParameterizationChanged();
}

double FaceGroup::AreaUV() const
{
    double areaUV = 0;
    for (auto fptr : fpVec) areaUV += std::abs(DistortionMetric::AreaUV(*fptr));
    return areaUV;
}

double FaceGroup::Area3D() const
{
    double area3D = 0;
    for (auto fptr : fpVec) area3D += DistortionMetric::Area3D(*fptr);
    return area3D;
}

double FaceGroup::BorderUV() const
{
    if (borderChanged) UpdateBorder();
    return borderUV;
}

vcg::Box2d FaceGroup::UVBox() const
{
    vcg::Box2d box;
    for (auto fptr : fpVec) {
        box.Add(fptr->WT(0).P());
        box.Add(fptr->WT(1).P());
        box.Add(fptr->WT(2).P());
    }
    return box;
}

void FaceGroup::ParameterizationChanged()
{
    borderChanged = true;
}

Mesh::FacePointer FaceGroup::Fp()
{
    assert(!fpVec.empty()); return fpVec[0];
}

std::size_t FaceGroup::FN() const
{
    return fpVec.size();
}

std::size_t FaceGroup::NumAdj() const
{
    return adj.size();
}

void FaceGroup::MapDistortion(DistortionMetric::Type type, ParameterizationGeometry geometry, double areaScale)
{
    minMappedFaceValue = std::numeric_limits<float>::max();
    maxMappedFaceValue = std::numeric_limits<float>::lowest();
    for (auto fptr : fpVec) {
        if (type == DistortionMetric::Area) {
            fptr->Q() = DistortionMetric::AreaDistortion(mesh, *fptr, areaScale, geometry);
        }
        else if (type == DistortionMetric::Angle) {
            fptr->Q() = DistortionMetric::AngleDistortion(mesh, *fptr, geometry);
        }
        else assert(0 && "FaceGroup::MapDistortion");
        minMappedFaceValue = std::min(minMappedFaceValue, fptr->Q());
        maxMappedFaceValue = std::max(maxMappedFaceValue, fptr->Q());
    }
}

void FaceGroup::UpdateBorder() const
{
    auto CCIDh = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(mesh, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(mesh, CCIDh));

    borderUV = 0.0;
    for (auto fptr : fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (face::IsBorder(*fptr, i) || CCIDh[fptr] != CCIDh[fptr->FFp(i)]) {
                borderUV += EdgeLengthUV(*fptr, i);
            }
        }
    }
    borderChanged = false;
}


// MeshGraph class implementation
// ==============================

MeshGraph::MeshGraph(Mesh& m)
    : mesh{m}
{
}

void MeshGraph::MapDistortion(DistortionMetric::Type type, ParameterizationGeometry geometry)
{
    double areaScale;
    DistortionMetric::ComputeAreaScale(mesh, areaScale, geometry);
    /// FIXME absolute distortion should be an option parameter
    if (geometry == ParameterizationGeometry::Texture) areaScale = 1;
    for (auto& c : charts) c.second->MapDistortion(type, geometry, areaScale);

    // Map distortion to color
    std::pair<float, float> range = DistortionRange();
    for (auto& c : charts) {
        for (auto fptr : c.second->fpVec) {
            float q = fptr->Q();
            //float q = fptr->Q() / std::max(std::abs(range.first), std::abs(range.second));
            //float q = fptr->Q() / std::abs(range.first);
            if (q < 0) {
                //float v = 1.0f - (q / range.first);
                float v = 1.0f + (q / (-range.first));
                fptr->C().Import(Color4f{1.0f, v, v, 1.0f});
            } else {
                //float v = 1.0f - (q / range.second);
                float v = 1.0f - (q / range.second);
                fptr->C().Import(Color4f{v, v, 1.0f, 1.0f});
            }
        }
    }
}

std::pair<float,float> MeshGraph::DistortionRange() const
{
    std::pair<float,float> range = std::make_pair(std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest());
    for (const auto& c : charts) {
        range.first = std::min(c.second->minMappedFaceValue, range.first);
        range.second = std::max(c.second->maxMappedFaceValue, range.second);
    }
    return range;
}

std::shared_ptr<FaceGroup> MeshGraph::GetChart(RegionID i)
{
    //assert(charts.find(i) != charts.end() && "Chart does not exist");
    auto e = charts.find(i);
    if (e != charts.end()) return e->second;
    else return nullptr;
}

std::shared_ptr<FaceGroup> MeshGraph::GetChart_Insert(RegionID i)
{
    if (charts.find(i) == charts.end()) charts.insert(std::make_pair(i, std::make_shared<FaceGroup>(mesh, i)));
    return charts[i];
}

std::size_t MeshGraph::Count() const
{
    return charts.size();
}

int MeshGraph::MergeCount() const
{
    int n = 0;
    for (const auto& c : charts) n += c.second->numMerges;
    return n;

}

double MeshGraph::Area3D() const
{
    double area3D = 0.0f;
    for (const auto& c : charts) area3D += c.second->Area3D();
    return area3D;
}

double MeshGraph::AreaUV() const
{
    double areaUV = 0.0f;
    for (const auto& c : charts) areaUV += c.second->AreaUV();
    return areaUV;
}

double MeshGraph::BorderUV(float *meshBorderLengthUV, float *seamLengthUV)
{
    double borderUV = 0.0;
    for (const auto& c : charts) borderUV += c.second->BorderUV();
    return borderUV;
}

// GraphManager class implementation
// =================================

GraphManager::GraphManager(std::shared_ptr<MeshGraph> gptr, std::unique_ptr<EdgeWeightFunction> &&wfct)
    : g{gptr},
      wfct{std::move(wfct)}
{
    tri::UpdateTopology<Mesh>::FaceFace(g->mesh);
    for (auto elem : g->charts) {
        ChartHandle c1 = elem.second;
        for (ChartHandle c2 : c1->adj) {
            AddEdge(c1, c2);
        }
    }
}


std::shared_ptr<MeshGraph> GraphManager::Graph()
{
    return g;
}

bool GraphManager::Valid(std::pair<GraphManager::Edge,double> weightedEdge)
{
    auto e = edges.find(weightedEdge.first);
    return e != edges.end() && e->second == weightedEdge.second;
}

const std::pair<GraphManager::Edge,double>& GraphManager::PeekNextEdge()
{
    //assert(HasNextEdge());
    return queue.top();
}

void GraphManager::RemoveNextEdge()
{
    //assert(HasNextEdge());
    queue.pop();
}

bool GraphManager::HasNextEdge()
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
                GraphManager::Edge e = we.first;
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
ChartHandle GraphManager::Collapse(const GraphManager::Edge& e)
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

bool GraphManager::AddEdge(ChartHandle c1, ChartHandle c2, bool replace)
{
     GraphManager::Edge e{c1, c2};

     if (replace) edges.erase(e);

     if (edges.find(e) == edges.end()) {
         auto weightedEdge = std::make_pair(e, (*wfct)(e));
         edges.insert(weightedEdge);
         queue.push(weightedEdge);
         return true;
     }
     else {
         return false;
     }
}

void GraphManager::Merge(ChartHandle c1, ChartHandle c2)
{
    assert(edges.find(GraphManager::Edge{c1,c2}) != edges.end());

    auto CCIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(g->mesh, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(g->mesh, CCIDh));

    // Update ID and add faces
    for (auto fp : c2->fpVec) {
        CCIDh[fp] = c1->id;
        c1->AddFace(fp);
    }

    // Update adjacencies - Note that obsolete edges remain in the priority queue
    c1->adj.erase(c2);
    for (auto cn : c2->adj) {
        GraphManager::Edge e2{c2, cn};
        edges.erase(e2);
        if (cn != c1) { // c1 is now (if it wasn't already) adjacent to cn
            cn->adj.erase(c2);
            cn->adj.insert(c1);
            c1->adj.insert(cn);

            // Manage queue and edge map state
            AddEdge(c1, cn, true); // replace edge
            /*
            GraphManager::Edge e1{c1, cn};
            edges.erase(e1); // delete the edge if it already exists so it can be reweighted
            auto weighted = std::make_pair(e1, wfct(e1));
            edges.insert(weighted);
            queue.push(weighted);
            */
        }
    }

    c1->numMerges += c2->numMerges + 1;

    g->charts.erase(c2->id);
}

void GraphManager::Split(const RegionID id, std::vector<ChartHandle>& splitCharts)
{
    assert(g->charts.count(id) == 1);

    auto CCIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(g->mesh, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(g->mesh, CCIDh));

    auto ICCh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(g->mesh, "InitialConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(g->mesh, ICCh));

    ChartHandle chart = g->GetChart(id);
    if (chart->numMerges == 0) {
        splitCharts.push_back(chart);
        return;
    }

    // detach chart from the graph
    g->charts.erase(chart->id);
    for (auto cn : chart->adj) {
        cn->adj.erase(chart); // remove reference in the neighbor
        edges.erase(GraphManager::Edge{chart, cn}); // remove edge from support list
    }

    // restore original ids
    for (auto fptr : chart->fpVec) {
        CCIDh[fptr] = ICCh[fptr];
    }

    std::set<ChartHandle> newCharts;
    // readd faces to the graph, creating the required nodes
    for (auto fptr : chart->fpVec) {
        RegionID regionId = CCIDh[fptr];
        auto ch = g->GetChart_Insert(regionId);
        newCharts.insert(ch);
        ch->AddFace(fptr);
        for (int i = 0; i < fptr->VN(); ++i) {
            RegionID adjId = CCIDh[fptr->FFp(i)];
            if (regionId != adjId) {
                // force symmetry
                (ch->adj).insert(g->GetChart_Insert(adjId));
                (g->GetChart_Insert(adjId)->adj).insert(ch);
            }
        }
    }

    // insert new edges in the support list
    for (auto c1 : newCharts) {
        splitCharts.push_back(c1);
        for (auto c2 : c1->adj) {
            AddEdge(c1, c2, true);
        }
    }
}

int GraphManager::CloseMacroRegions(std::size_t minRegionSize)
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
