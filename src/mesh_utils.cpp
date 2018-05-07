#include "mesh_utils.h"

#include <vcg/complex/complex.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/complex/algorithms/geodesic.h>

#include <vector>
#include <algorithm>

using namespace vcg;

struct PosNode {
    PosF pos;
    double distance;

    PosNode() = default;
    PosNode(PosF p, double d) : pos{p}, distance{d} {}

    bool operator<(const PosNode& other) { return distance < other.distance; }
};

// assumption: the original texture seams are marked as crease edges
double ComputeDistanceFromBorderOnSeams(Mesh& m)
{
    tri::UpdateTopology<Mesh>::FaceFace(m);
    tri::UpdateFlags<Mesh>::VertexBorderFromFaceAdj(m);
    tri::UpdateQuality<Mesh>::VertexConstant(m, INFINITY);

    // initialize the frontier as a sequence of PosNode objects pointing to boundary vertices along texture seams
    std::vector<PosNode> probes;
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            if (f.IsCrease(i)) {
                PosF p{&f, i};
                if (p.V()->IsB() && p.V()->Q() == INFINITY) {
                    probes.push_back(PosNode{p, 0});
                    p.V()->Q() = 0;
                }
                p.FlipV();
                if (p.V()->IsB() && p.V()->Q() == INFINITY) {
                    p.V()->Q() = 0;
                    probes.push_back(PosNode{p, 0});
                }
            }
        }
    }

    tri::EuclideanDistance<Mesh> dist;
    std::make_heap(probes.begin(), probes.end());
    while (!probes.empty()) {
        std::pop_heap(probes.begin(), probes.end());
        PosNode node = probes.back();
        probes.pop_back();
        if (node.distance == node.pos.V()->Q()) {
            std::vector<PosF> fan = GetCreasePosFan(node.pos);
            for (auto& creasePos : fan) {
                double d = node.pos.V()->Q() + dist(node.pos.V(), creasePos.V());
                if (d < creasePos.V()->Q()) {
                    creasePos.V()->Q() = d;
                    probes.push_back(PosNode{creasePos, d});
                    std::push_heap(probes.begin(), probes.end());
                }
            }
        } // otherwise the entry is obsolete
    }

    double maxDist = 0;
    for (auto& v : m.vert) {
        if (v.Q() < 0) assert(0);
        if (v.Q() < INFINITY && v.Q() > maxDist) maxDist = v.Q();
    }

    return maxDist;
}


void MarkInitialSeamsAsCreases(Mesh& m, SimpleTempData<Mesh::FaceContainer, RegionID>& initialId)
{
    tri::UpdateFlags<Mesh>::FaceClearCreases(m);
    for (auto& f : m.face) {
        f.C() = vcg::Color4b::DarkGray;
    }
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            if (!f.IsB(i)) {
                auto& ff = *(f.FFp(i));
                int ii = f.FFi(i);
                auto p = std::make_pair(tri::Index(m, f), tri::Index(m, ff));
                if (!f.IsCrease(i) && initialId[f] != initialId[ff]) {
                    f.SetCrease(i);
                    ff.SetCrease(ii);
                    f.C() = vcg::Color4b::Blue;
                    ff.C() = vcg::Color4b::Blue;
                }
            }
        }
    }
}

void ClearCreaseLoops(Mesh& m)
{
    for (auto& f : m.face) {
        for (int i = 0; i < 3; i++) {
            if (f.IsCrease(i)) {

            }
        }
    }
}

std::vector<PosF> GetCreasePosFan(PosF& startPos)
{
    std::vector<PosF> posVec;
    PosF p = startPos;
    do {
        if (p.IsCrease()) {
            PosF pp = p;
            pp.FlipV();
            posVec.push_back(pp);
        }
        p.FlipF();
        p.FlipE();
    } while (p != startPos);
    return posVec;
}

void ClearCreasesAlongPath(Mesh& m, std::vector<PosF>& path)
{
    for (auto p : path) {
        Mesh::FacePointer startF = p.F();
        p.F()->ClearCrease(p.E());
        p.FlipF();
        p.F()->ClearCrease(p.E());
        assert(p.F() == startF && "Mesh is non-manifold along path");
        assert(!p.IsBorder());
        assert(!p.V()->IsB());
    }
}

void MarkPathAsNonFaux(Mesh& m, std::vector<PosF>& path)
{
    tri::UpdateFlags<Mesh>::FaceSetF(m);
    for (auto p : path) {
        Mesh::FacePointer startF = p.F();
        p.F()->ClearF(p.E());
        p.FlipF();
        p.F()->ClearF(p.E());
        assert(p.F() == startF && "Mesh is non-manifold along path");
    }
}

// assumes initial texture seams are marked as creases
// the computed path is returned as a sequence of pos objects in the 'path' vector
// returns true if the path reaches the boundary, false if it loops (it can happen if islands are merged)
bool ComputePathToBoundaryAlongSeams(Mesh& m, const PosF& pos, std::vector<PosF>& path)
{
    tri::UpdateFlags<Mesh>::VertexClearV(m);
    path.clear();


    // each 'node' is a pair <Pos, d>
    // the Pos references a vertex and edge on the path
    std::unordered_map<Mesh::VertexPointer, PosNode> P; // predecessors

    // the assumption is that the initial pos points to the source node and edge
    // from where the edge starts

    // the predecessor of the source is the source itself
    P.insert(std::make_pair(pos.V(), PosNode(pos, 0)));

    std::deque<PosNode> Q;
    Q.push_back(PosNode{pos, 0});

    tri::EuclideanDistance<Mesh> dist;

    bool foundBorder = false;
    PosF marker;
    while (!Q.empty()) {
        std::sort(Q.begin(), Q.end());
        PosNode node = Q.front();
        Q.pop_front();

        PosF& p = node.pos;

        assert(P.count(p.V()) != 0);
        assert(p.F()->IsCrease(p.E()));

        p.V()->SetV();
        marker = p;

        if (p.V()->IsB()) {
            foundBorder = true;
            break;
        } else if (p.V() == pos.V()) {
            // in this case we looped back to the beginning, it means that the seam does not reach the mesh boundary
            foundBorder = false;
            break;
        }

        std::vector<PosF> creaseFan = GetCreasePosFan(p);
        for (auto& pcf : creaseFan) {
            double d = node.distance + dist(p.V(), pcf.V());
            auto pred = P.find(pcf.V());
            if (pred == P.end() || pred->second.distance > d) {
                PosNode nn{pcf, d};
                P[pcf.V()] = node; // node becomes the predecessor of the new node nn
                if (pred != P.end()) { // before inserting the node in the queue, remove the old entry if it exists
                    Q.erase(std::remove_if(Q.begin(), Q.end(), [&p](const PosNode& elem) { return elem.pos.V() == p.V(); }),
                            Q.end());
                }
                Q.push_back(nn);
            }
        }
    }

    // walk back from the frontier to the source, building the path (a sequence of Pos objects)

    while (marker.V() != pos.V()) {
        std::cout << "traversing face " << tri::Index(m, marker.F()) << " vertex " << marker.VInd() << std::endl;
        path.push_back(marker);
        marker = P[marker.V()].pos;
    }

}

// assumes initial texture seams are marked as creases
void MarkShortestSeamToBorderAsNonFaux(Mesh& m, const PosF& pos)
{
    tri::UpdateFlags<Mesh>::VertexClearV(m);


    // each 'node' is a pair <Pos, d>
    // the Pos references a vertex and edge on the path
    std::unordered_map<Mesh::VertexPointer, PosNode> P; // predecessors

    // the assumption is that the initial pos points to the source node and edge
    // from where the edge starts

    // the predecessor of the source is the source itself
    P.insert(std::make_pair(pos.V(), PosNode(pos, 0)));

    std::deque<PosNode> Q;
    Q.push_back(PosNode{pos, 0});

    tri::EuclideanDistance<Mesh> dist;

    PosF cutter;
    while (!Q.empty()) {
        std::sort(Q.begin(), Q.end());
        PosNode node = Q.front();
        Q.pop_front();

        PosF& p = node.pos;

        assert(P.count(p.V()) != 0);
        assert(p.F()->IsCrease(p.E()));

        p.V()->SetV();
        if (p.V()->IsB()) {
            cutter = p;
            break;
        }

        std::vector<PosF> creaseFan = GetCreasePosFan(p);
        for (auto& pcf : creaseFan) {
            double d = node.distance + dist(p.V(), pcf.V());
            auto pred = P.find(pcf.V());
            if (pred == P.end() || pred->second.distance > d) {
                PosNode nn{pcf, d};
                P[pcf.V()] = node; // node becomes the predecessor of the new node nn
                if (pred != P.end()) { // before inserting the node in the queue, remove the old entry if it exists
                    Q.erase(std::remove_if(Q.begin(), Q.end(), [&p](const PosNode& elem) { return elem.pos.V() == p.V(); }),
                            Q.end());
                }
                Q.push_back(nn);
            }
        }
    }

    // walk back from the frontier to the source, marking the edges as non faux
    // edges are marked on each face otherwise the cutting function does not work
    tri::UpdateFlags<Mesh>::FaceSetF(m);


    while (cutter.V() != pos.V()) {
       // std::cout << "traversing face " << tri::Index(m, cutter.F()) << " vertex " << cutter.VInd() << std::endl;
        assert(P.count(cutter.V()) == 1);
        Mesh::FacePointer fp = cutter.F();
        cutter.F()->ClearF(cutter.E());
        cutter.FlipF();
        cutter.F()->ClearF(cutter.E());
        cutter.FlipF();
        assert(cutter.F() == fp);
        PosNode& pn = P[cutter.V()];
        assert(!pn.pos.IsBorder());
        cutter = pn.pos;
    }

}
