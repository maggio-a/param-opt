#include "mesh_utils.h"
#include "mesh_attribute.h"

#include <vcg/complex/complex.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/complex/algorithms/geodesic.h>

#include <vector>
#include <algorithm>

using namespace vcg;

void MarkInitialSeamsAsFaux(Mesh& shell, Mesh& baseMesh)
{
    tri::UpdateTopology<Mesh>::FaceFace(shell);
    tri::UpdateFlags<Mesh>::FaceBorderFromFF(shell);
    tri::UpdateFlags<Mesh>::FaceClearF(shell);
    assert(HasFaceIndexAttribute(shell));
    auto ia = GetFaceIndexAttribute(shell);
    for (auto& sf : shell.face) {
        if (sf.holeFilling == false) {
            for (int i = 0; i < 3; ++i) {
                auto& f = baseMesh.face[ia[sf]];
                if (f.IsF(i))
                    sf.SetF(i);
            }
        }
    }
}

std::vector<PosF> GetFauxPosFan(PosF& startPos)
{
    std::vector<PosF> posVec;
    PosF p = startPos;
    do {
        if (p.IsFaux()) {
            PosF pp = p;
            pp.FlipV();
            posVec.push_back(pp);
        }
        p.FlipF();
        p.FlipE();
    } while (p != startPos);
    return posVec;
}

/* Function object to compute the average length across an edge, according
 * to the target shape of the shell faces. This allows to compute meaningful
 * seam path lengths on shell objects */
struct FeatureBasedEdgeLength {
    Mesh::PerFaceAttributeHandle<CoordStorage> targetShape;

    double operator()(PosF p) {
        CoordStorage cs = targetShape[p.F()];
        double l1 = (cs.P[p.E()] - cs.P[(p.E()+1)%3]).Norm();
        p.FlipF();
        double l2 = (cs.P[p.E()] - cs.P[(p.E()+1)%3]).Norm();
        return (l1 + l2) / 2.0;
    };

};

double ComputeDistanceFromBorderOnSeams(Mesh& m)
{
    tri::UpdateTopology<Mesh>::FaceFace(m);
    tri::UpdateFlags<Mesh>::VertexBorderFromFaceAdj(m);
    tri::UpdateQuality<Mesh>::VertexConstant(m, INFINITY);

    // initialize the frontier as a sequence of PosNode objects pointing to boundary vertices along texture seams
    std::vector<PosNode> probes;
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            if (f.IsF(i)) {
                PosF p{&f, i};
                if (p.V()->IsB() && p.V()->Q() == INFINITY) {
                    probes.push_back(PosNode{p, 0});
                    p.V()->Q() = 0;
                }
                p.FlipV();
                if (p.V()->IsB() && p.V()->Q() == INFINITY) {
                    probes.push_back(PosNode{p, 0});
                    p.V()->Q() = 0;
                }
            }
        }
    }

    assert(HasTargetShapeAttribute(m));
    auto targetShape = GetTargetShapeAttribute(m);
    FeatureBasedEdgeLength dist{targetShape};
    auto posNodeComp = [] (const PosNode& a, const PosNode& b) { return a.distance > b.distance; };
    //std::make_heap(probes.begin(), probes.end(), posNodeComp);
    while (!probes.empty()) {
        std::pop_heap(probes.begin(), probes.end(), posNodeComp);
        PosNode node = probes.back();
        probes.pop_back();
        if (node.distance == node.pos.V()->Q()) {
            std::vector<PosF> fan = GetFauxPosFan(node.pos);
            for (auto& fauxPos : fan) {
                double d = node.pos.V()->Q() + dist(fauxPos);
                assert(d > node.pos.V()->Q());
                if (d < fauxPos.V()->Q()) {
                    assert(fauxPos.V()->IsB() == false);
                    fauxPos.V()->Q() = d;
                    probes.push_back(PosNode{fauxPos, d});
                    std::push_heap(probes.begin(), probes.end(), posNodeComp);
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

void SelectShortestSeamPathToBoundary(Mesh& m, const PosF& pos)
{
    assert(pos.IsFaux());

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
        assert(p.F()->IsF(p.E()));

        p.V()->SetV();
        if (p.V()->IsB()) {
            cutter = p;
            break;
        }

        std::vector<PosF> fauxFan = GetFauxPosFan(p);
        for (auto& pcf : fauxFan) {
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
    tri::UpdateFlags<Mesh>::FaceClearFaceEdgeS(m);

    while (cutter.V() != pos.V()) {
        assert(P.count(cutter.V()) == 1);
        Mesh::FacePointer fp = cutter.F();
        cutter.F()->SetFaceEdgeS(cutter.E());
        cutter.FlipF();
        cutter.F()->SetFaceEdgeS(cutter.E());
        cutter.FlipF();
        assert(cutter.F() == fp && "Mesh is not edge manifold along path");
        PosNode& pn = P[cutter.V()];
        assert(!pn.pos.IsBorder());
        cutter = pn.pos;
    }
}

void SelectShortestSeamPathToPeak(Mesh& m, const PosF& pos)
{
    assert(pos.IsFaux());
    tri::UpdateFlags<Mesh>::VertexClearV(m);
    PosF curPos = pos;
    if (curPos.V()->Q() < curPos.VFlip()->Q())
        curPos.FlipV();
    curPos.VFlip()->SetV();
    double curDistance = pos.V()->Q();
    auto posQualityComparator = [](const PosF& p1, const PosF& p2) {
        return p1.V()->Q() < p2.V()->Q();
    };
    while (true) {
        Mesh::FacePointer fp = curPos.F();
        assert(curPos.VFlip()->IsV());
        curPos.F()->SetFaceEdgeS(curPos.E());
        curPos.FlipF();
        curPos.F()->SetFaceEdgeS(curPos.E());
        curPos.FlipF();
        assert(curPos.F() == fp && "Mesh is not edge manifold along path");
        curPos.V()->SetV();
        auto pred = [&](const PosF& p){
            return p.V()->IsV() || (p.V()->Q() < curDistance);
        };
        std::vector<PosF> fan = GetFauxPosFan(curPos);
        fan.erase(std::remove_if(fan.begin(), fan.end(), pred), fan.end());
        if (fan.empty())
            break;
        else {
            //curPos = *(std::min_element(fan.begin(), fan.end(), posQualityComparator));
            curPos = *(std::max_element(fan.begin(), fan.end(), posQualityComparator));
            assert(curDistance <= curPos.V()->Q());
            curDistance = curPos.V()->Q();
        }
    }
}

void CleanupShell(Mesh& shell)
{
    tri::UpdateTopology<Mesh>::FaceFace(shell);
    tri::UpdateFlags<Mesh>::VertexBorderFromFaceAdj(shell);
    for (auto& sf : shell.face) {
        for (int i = 0; i < 3; i++) {
            if (face::IsBorder(sf, i)) sf.ClearF(i);
        }
    }
    tri::UpdateFlags<Mesh>::FaceClearV(shell);
    for (auto& sf : shell.face) {
        if (sf.holeFilling && !sf.IsV()) {
            bool boundary = false;
            for (int i = 0; i < 3; ++i) {
                if (face::IsBorder(sf, i)) boundary = true;
            }
            if (boundary) {
                // vist the holeFilling region
                std::stack<Mesh::FacePointer> s;
                s.push(&sf);
                while (!s.empty()) {
                    Mesh::FacePointer fp = s.top();
                    s.pop();
                    fp->SetV();
                    for (int i = 0; i < 3; ++i) {
                        if (!fp->FFp(i)->IsV() && fp->FFp(i)->holeFilling) {
                            s.push(fp->FFp(i));
                        }
                    }
                }
            }
        }
    }
    for (auto& sf : shell.face) {
        if (sf.IsV()) tri::Allocator<Mesh>::DeleteFace(shell, sf);
    }
    int removed = tri::Clean<Mesh>::RemoveUnreferencedVertex(shell);
    if (removed > 0) {
        std::cout << removed << " unreferenced vertices removed after cutting" << std::endl;
        tri::Allocator<Mesh>::CompactVertexVector(shell);
    }
    tri::Allocator<Mesh>::CompactEveryVector(shell);
}

