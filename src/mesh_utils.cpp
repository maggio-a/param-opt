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
    ensure_condition(HasFaceIndexAttribute(shell));
    auto ia = GetFaceIndexAttribute(shell);
    for (auto& sf : shell.face) {
        if (sf.IsMesh()) {
            for (int i = 0; i < 3; ++i) {
                auto& f = baseMesh.face[ia[sf]];
                if (f.IsF(i) && !face::IsBorder(sf, i)) {
                    sf.SetF(i);
                    sf.FFp(i)->SetF(sf.FFi(i));
                }
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
        if (p.IsBorder()) break;
    } while (p != startPos);

    if (p.IsBorder() && p != startPos) {
        p = startPos;
        p.FlipE();
        while (true) {
            if (p.IsFaux()) {
                PosF pp = p;
                pp.FlipV();
                posVec.push_back(pp);
            }
            if (p.IsBorder()) break;
            p.FlipF();
            p.FlipE();
            ensure_condition(p != startPos);
        }
    }

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

struct VertexPosEdgeLength {
    double operator()(PosF p) {
        return (p.V()->P() - p.VFlip()->P()).Norm();
    }
};

void ComputeDistanceFromBorderOnSeams(Mesh& m)
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

    ensure_condition(HasTargetShapeAttribute(m));
    auto targetShape = GetTargetShapeAttribute(m);
    //FeatureBasedEdgeLength dist{targetShape};

    VertexPosEdgeLength dist;


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
                ensure_condition(d > node.pos.V()->Q());
                if (d < fauxPos.V()->Q()) {
                    ensure_condition(fauxPos.V()->IsB() == false);
                    fauxPos.V()->Q() = d;
                    probes.push_back(PosNode{fauxPos, d});
                    std::push_heap(probes.begin(), probes.end(), posNodeComp);
                }
            }
        } // otherwise the entry is obsolete
    }
}

PosF SelectShortestSeamPathToBoundary(Mesh& m, const PosF& pos)
{
    ensure_condition(pos.IsFaux());
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
        ensure_condition(curPos.VFlip()->IsV());
        curPos.F()->SetFaceEdgeS(curPos.E());
        curPos.FlipF();
        curPos.F()->SetFaceEdgeS(curPos.E());
        curPos.FlipF();
        ensure_condition(curPos.F() == fp && "Mesh is not edge manifold along path");
        curPos.V()->SetV();
        if (curPos.V()->IsB())
            break;
        else {
            std::vector<PosF> fan = GetFauxPosFan(curPos);
            curPos = *(std::min_element(fan.begin(), fan.end(), posQualityComparator));
            ensure_condition(curDistance >= curPos.V()->Q());
            ensure_condition(!curPos.V()->IsV());
            curDistance = curPos.V()->Q();
        }
    }
    return curPos;
}

void SelectShortestSeamPathToPeak(Mesh& m, const PosF& pos)
{
    ensure_condition(pos.IsFaux());
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
        ensure_condition(curPos.VFlip()->IsV());
        curPos.F()->SetFaceEdgeS(curPos.E());
        curPos.FlipF();
        curPos.F()->SetFaceEdgeS(curPos.E());
        curPos.FlipF();
        ensure_condition(curPos.F() == fp && "Mesh is not edge manifold along path");
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
            ensure_condition(curDistance <= curPos.V()->Q());
            curDistance = curPos.V()->Q();
        }
    }
}

/* Start from a boundary pos and walk on the selected path. If at any point, a
 * vertex path adjacent to a hole-filling face is found, from that vertex onward
 * clear the path. This allows to compute paths that join chart boundaries
 * without accidentally cutting away small shell faces that then remain isolated
 * in the parameterization */
bool RectifyCut(Mesh& shell, PosF boundaryPos)
{
    ensure_condition(boundaryPos.V()->IsB());
    tri::UpdateFlags<Mesh>::VertexClearV(shell);
    PosF p = boundaryPos;
    bool reachedFillArea = false;
    while (true) {
        p.V()->SetV();

        if (reachedFillArea) {
            p.F()->ClearFaceEdgeS(p.E());
            p.FlipF();
            p.F()->ClearFaceEdgeS(p.E());
            p.FlipF();
        }

        if (!p.IsBorder() && !reachedFillArea) {
            PosF pp = p;
            do {
                if (pp.F()->IsHoleFilling())
                    reachedFillArea = true;
                pp.FlipF();
                pp.FlipE();
            } while (p != pp && !reachedFillArea);
        }

        std::vector<PosF> fan = GetFauxPosFan(p);
        bool moved = false;
        for (auto& fanPos : fan) {
            if (!fanPos.V()->IsV()) {
                if (fanPos.F()->IsFaceEdgeS(fanPos.E())) {
                    ensure_condition(!moved && "Cut path has branches");
                    moved = true;
                    p = fanPos;
                }
            }
        }
        if (!moved) break;
    }
    return reachedFillArea;
}

void CleanupShell(Mesh& shell)
{
    // Select vertices that were incident on the cut path
    tri::UpdateFlags<Mesh>::VertexClearS(shell);
    for (auto& sf : shell.face) {
        for (int i = 0; i < 3; ++i) {
            if (sf.IsFaceEdgeS(i)) {
                sf.V0(i)->SetS();
                sf.V1(i)->SetS();
            }
        }
    }
    // If a fill area touches the cut even with one single vertex, it must be removed
    tri::UpdateFlags<Mesh>::FaceClearV(shell);
    for (auto& sf : shell.face) {
        if (sf.IsHoleFilling() && !sf.IsV()) {
            bool fillAreaReached = false;
            for (int i = 0; i < 3; ++i) {
                if (sf.V(i)->IsS()) fillAreaReached = true;
            }
            if (fillAreaReached) {
                // vist the holeFilling region
                std::stack<Mesh::FacePointer> s;
                s.push(&sf);
                while (!s.empty()) {
                    Mesh::FacePointer fp = s.top();
                    s.pop();
                    fp->SetV();
                    for (int i = 0; i < 3; ++i) {
                        if (!fp->FFp(i)->IsV() && fp->FFp(i)->IsHoleFilling()) {
                            s.push(fp->FFp(i));
                        }
                    }
                }
            }
        }
    }

    bool deleted = false;
    for (auto& sf : shell.face) {
        if (sf.IsV()) {
            tri::Allocator<Mesh>::DeleteFace(shell, sf);
            deleted = true;
        }
    }
    tri::Clean<Mesh>::RemoveUnreferencedVertex(shell);
    tri::Allocator<Mesh>::CompactEveryVector(shell);

    tri::UpdateTopology<Mesh>::FaceFace(shell);
    for (auto& sf : shell.face) {
        for (int i = 0; i < 3; i++) {
            if (face::IsBorder(sf, i)) sf.ClearF(i);
        }
    }

    // It can happen that after deleting the hole-filling faces the mesh is
    // no longer vertex manifold because the new boundary touched only a
    // vertex of a hole-filling face
    if (deleted) {
        tri::Clean<Mesh>::SplitNonManifoldVertex(shell, 0.15);
    }
}

void CopyShell(Mesh& shell, Mesh& out)
{
    out.Clear();
    out.ClearAttributes();

    ensure_condition(HasBoundaryInfoAttribute(shell));
    ensure_condition(HasTargetShapeAttribute(shell));
    ensure_condition(HasFaceIndexAttribute(shell));
    GetBoundaryInfoAttribute(out);
    GetTargetShapeAttribute(out);
    GetFaceIndexAttribute(out);

    tri::Append<Mesh,Mesh>::MeshCopy(out, shell, false, true);
}
