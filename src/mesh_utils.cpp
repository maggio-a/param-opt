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
    tri::UpdateFlags<Mesh>::FaceClearF(shell);
    assert(HasFaceIndexAttribute(shell));
    auto ia = GetFaceIndexAttribute(shell);
    for (auto& sf : shell.face) {
        if (sf.holeFilling == false) {
            for (int i = 0; i < 3; ++i) {
                auto& f = baseMesh.face[ia[sf]];
                if (f.IsF(i) && !face::IsBorder(sf, i))
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
    tri::UpdateFlags<Mesh>::VertexBorderFromFaceAdj(m);
    tri::UpdateFlags<Mesh>::VertexClearV(m);
    PosF curPos = pos;
    if (curPos.V()->Q() > curPos.VFlip()->Q())
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
        if (curPos.V()->IsB())
            break;
        else {
            std::vector<PosF> fan = GetFauxPosFan(curPos);
            curPos = *(std::min_element(fan.begin(), fan.end(), posQualityComparator));
            assert(curDistance >= curPos.V()->Q());
            assert(!curPos.V()->IsS());
            curDistance = curPos.V()->Q();
        }
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
    tri::Clean<Mesh>::RemoveUnreferencedVertex(shell);
    tri::Allocator<Mesh>::CompactEveryVector(shell);
}

void CopyShell(Mesh& shell, Mesh& out)
{
    out.Clear();
    out.ClearAttributes();

    assert(HasBoundaryInfoAttribute(shell));
    assert(HasTargetShapeAttribute(shell));
    assert(HasFaceIndexAttribute(shell));
    GetBoundaryInfoAttribute(out);
    GetTargetShapeAttribute(out);
    GetFaceIndexAttribute(out);

    tri::Append<Mesh,Mesh>::MeshCopy(out, shell, false, true);
}
