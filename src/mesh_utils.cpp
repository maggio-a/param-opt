#include "mesh_utils.h"
#include "math_utils.h"
#include "mesh_attribute.h"
#include "logging.h"
#include "timer.h"

#include <vcg/complex/complex.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/complex/algorithms/geodesic.h>

#include <vcg/complex/algorithms/isotropic_remeshing.h>

#include <vector>
#include <algorithm>


using namespace vcg;


// assumes topology is updated (FaceFace)
void ComputeBoundaryInfo(Mesh& m)
{
    BoundaryInfo& info = GetBoundaryInfoAttribute(m)();
    info.Clear();
    tri::UpdateFlags<Mesh>::FaceClearV(m);
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            if (!f.IsV() && face::IsBorder(f, i)) {
                double totalBorderLength = 0;
                std::vector<std::size_t> borderFaces;
                std::vector<int> vi;

                face::Pos<Mesh::FaceType> p(&f, i);
                face::Pos<Mesh::FaceType> startPos = p;
                ensure_condition(p.IsBorder());
                do {
                    ensure_condition(p.IsManifold());
                    p.F()->SetV();
                    borderFaces.push_back(tri::Index(m, p.F()));
                    vi.push_back(p.VInd());
                    totalBorderLength += EdgeLength(*p.F(), p.VInd());
                    p.NextB();
                } while (p != startPos);
                info.vBoundaryLength.push_back(totalBorderLength);
                info.vBoundarySize.push_back(borderFaces.size());
                info.vBoundaryFaces.push_back(borderFaces);
                info.vVi.push_back(vi);
            }
        }
    }

    LOG_DEBUG << "Mesh has " << info.N() << " boundaries";
}

void CloseMeshHoles(Mesh& shell)
{
    static double t0 = 0;
    static double t1 = 0;
    static double t2 = 0;
    static double t3 = 0;
    static double t4 = 0;
    static double t5 = 0;

    Timer t;

    int startFN = shell.FN();

    // Get border info
    ComputeBoundaryInfo(shell);
    BoundaryInfo& info = GetBoundaryInfoAttribute(shell)();

    t0 += t.TimeElapsed();

    //vcg::tri::io::ExporterPLY<Mesh>::Save(shell, "original.ply", tri::io::Mask::IOM_FACECOLOR);
    // Leave only the longest boundary
    tri::UpdateFlags<Mesh>::FaceClearS(shell);
    ensure_condition(info.vBoundaryFaces.size() > 0 && "Mesh has no boundaries");
    if (info.vBoundaryFaces.size() > 1) {
        std::size_t k = info.LongestBoundary();
        // select all the boundary faces
        for (std::size_t i = 0; i < info.vBoundaryFaces.size(); ++i) {
            if (i == k) continue;
            for (auto j : info.vBoundaryFaces[i]) {
                ensure_condition(face::IsBorder(shell.face[j], 0) || face::IsBorder(shell.face[j], 1) || face::IsBorder(shell.face[j], 2));
                shell.face[j].SetS();
            }
        }
        tri::Hole<Mesh>::EarCuttingFill<tri::MinimumWeightEar<Mesh>>(shell, shell.FN(), true);
    }

    t1 += t.TimeElapsed();

    // Compute border info
    double length = 0;
    int count = 0;
    for (std::size_t i = 0; i < info.vBoundaryFaces.size(); ++i) {
        if (i == info.LongestBoundary())
            continue;
        length += info.vBoundaryLength[i];
        count += (int) info.vBoundarySize[i];
    }

    // select the hole-filling faces that need to be remeshed
    tri::UpdateFlags<Mesh>::FaceClearS(shell);
    for (auto& sf : shell.face) {
        if (sf.IsHoleFilling())
            sf.SetS();
    }

    t2 += t.TimeElapsed();

    // Remesh filled hole
    //ColorFace(shell);
    IsotropicRemeshing<Mesh>::Params params;
    //params.SetTargetLen(2.0*(totalBorderLen / totalBorderFaces));
    //params.SetTargetLen(length / count);
    params.SetTargetLen(length / (count));
    params.SetFeatureAngleDeg(30);
    params.selectedOnly = true;
    //params.splitFlag = false;
    //params.swapFlag = false;
    //params.collapseFlag = false;
    params.smoothFlag = false;
    params.projectFlag = false;
    params.iter = 3;
    int iter = 0;
    do {
        IsotropicRemeshing<Mesh>::Do(shell, params);
        LOG_DEBUG << "Remeshing: " << params.stat.collapseNum << " " << params.stat.flipNum << " " << params.stat.splitNum;
        iter += params.iter;
    } while (params.stat.collapseNum + params.stat.flipNum + params.stat.splitNum > 0 && iter < 3);

    t3 += t.TimeElapsed();

    tri::Allocator<Mesh>::CompactEveryVector(shell);

    auto ia = GetFaceIndexAttribute(shell);
    tri::UpdateFlags<Mesh>::FaceClearS(shell);
    for (auto& f: shell.face) {
        if (int(tri::Index(shell, f)) >= startFN) {
            f.SetHoleFilling();
            ia[f] = -1;
        }
    }

    tri::UpdateTopology<Mesh>::FaceFace(shell);
    tri::UpdateTopology<Mesh>::VertexFace(shell);

    t4 += t.TimeElapsed();

    LOG_DEBUG << "\t\t\t[TIMING CloseMeshHoles()] t0 = " << t0;
    LOG_DEBUG << "\t\t\t[TIMING CloseMeshHoles()] t1 = " << t1;
    LOG_DEBUG << "\t\t\t[TIMING CloseMeshHoles()] t2 = " << t2;
    LOG_DEBUG << "\t\t\t[TIMING CloseMeshHoles()] t3 = " << t3;
    LOG_DEBUG << "\t\t\t[TIMING CloseMeshHoles()] t4 = " << t4;
}

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
    tri::UpdateQuality<Mesh>::VertexConstant(m, Infinity());

    // initialize the frontier as a sequence of PosNode objects pointing to boundary vertices along texture seams
    std::vector<PosNode> probes;
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            if (f.IsF(i)) {
                PosF p{&f, i};
                if (p.V()->IsB() && p.V()->Q() == Infinity()) {
                    probes.push_back(PosNode{p, 0});
                    p.V()->Q() = 0;
                }
                p.FlipV();
                if (p.V()->IsB() && p.V()->Q() == Infinity()) {
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
void RectifyCut(Mesh& shell, PosF boundaryPos)
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
        int nv;
        while ((nv = tri::Clean<Mesh>::SplitNonManifoldVertex(shell, 0.15)) > 0)
               ;
        tri::UpdateTopology<Mesh>::FaceFace(shell);
        tri::UpdateTopology<Mesh>::VertexFace(shell);
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

#include "uv.h"

std::vector<PosF> ExtractFauxTrack(Mesh& m, PosF pos)
{
    LOG_INFO << "Extracting track";
    std::vector<PosF> track;

    PosF p = pos;

    {
        int ind = p.VInd();
        p.FlipV();
        if (p.VInd() != ((ind + 1)%3))
            p.FlipV();
    }


    PosF start = p;
    do {
        ensure_condition(p.IsFaux());
        track.push_back(p);
        p.FlipE();
        while (!p.IsFaux()) {
            p.FlipF();
            p.FlipE();
        }
        p.FlipV();
    } while (p != start);


    /*
    PosF start = p;
    do {
        ensure_condition(p.IsFaux());
        track.push_back(p);
        do {
           p.NextE();
        } while (!p.IsFaux());
        p.FlipF();
        p.FlipV();
    } while (p != start);
    */

    LOG_INFO << "Extracted track of length " << track.size();

    return track;
}

#include <unordered_map>

#include <wrap/io_trimesh/export.h>

std::pair<VPU,FPU> PadFauxTrack(Mesh& m, std::vector<PosF>& track)
{
    ensure_condition(track.size() > 0);


    MyMesh em;
    for (auto pos : track)
        tri::Allocator<MyMesh>::AddEdge(em, pos.V()->P(), pos.VFlip()->P());

    tri::Clean<MyMesh>::RemoveDuplicateEdge(em);
    tri::Clean<MyMesh>::RemoveDuplicateVertex(em);
    tri::Allocator<MyMesh>::CompactEveryVector(em);

    tri::io::Exporter<MyMesh>::Save(em, "track.obj", tri::io::Mask::IOM_ALL);


    int nv = track.size();
    int nf = 2 * track.size();

    LOG_INFO << "PadFauxTrack(): Adding " << nv << " vertices and " << nf << " faces...";

    VPU vpu;
    FPU fpu;

    auto vi = tri::Allocator<Mesh>::AddVertices(m, nv, vpu);
    auto fi = tri::Allocator<Mesh>::AddFaces(m, nf, fpu);

    for (auto& pos : track) {
        vpu.Update(pos.v);
        fpu.Update(pos.f);
    }

    /* I do the padding in two passes,
     * First, I iterate over the track and duplicate the vertices while also
     *  building an index of each duplicate vertex pair. During this pass I also
     *  update the vertex references in the 'inner' fan, effectively detaching the
     *  chart outlined by the pos track from the rest of the mesh primitives.
     * Then, I iterate over the pos objects once again, creating the necessary faces
     * and (TODO) updating the face-face topology information
     */

    // first pass
    std::unordered_map<Mesh::VertexPointer, Mesh::VertexPointer> vmap;
    std::unordered_map<Mesh::VertexPointer, Mesh::VertexPointer> vmap_inv;

    for (auto& pos : track) {
        // retrieve the 'inner fan'
        std::vector<std::pair<Mesh::FacePointer, int>> innerfan;
        PosF p = pos;
        while (true) {
            innerfan.push_back(std::make_pair(p.F(), p.VInd()));
            p.FlipE();
            if (p.IsFaux())
                break;
            else
                p.FlipF();
        }

        // duplicate the vertex
        Mesh::VertexPointer vp = pos.V();
        //ensure_condition(vmap.count(vp) == 0); FIXME TODO this is triggered by non vertex-manifold charts
        if (vmap.count(vp) == 0) {
            vi->ImportData(*vp);
            vmap[vp] = &*vi;
            vi++;
        }
        Mesh::VertexPointer vpdup = vmap[vp];

        // critical, update the vertex pointer of the pos
        pos.V() = vpdup;

        vmap_inv[vpdup] = vp;

        // detach the faceface topology
        if (!face::IsBorder(*(pos.F()), pos.E()))
            vcg::face::FFDetach(*(pos.F()), pos.E());

        for (auto fanelement : innerfan) {
            Mesh::FacePointer fp = fanelement.first;
            int vind = fanelement.second;
            fp->V(vind) = vpdup;
        }
    }

    // second pass
    for (auto pos : track) {
        Mesh::VertexPointer vp1dup = pos.VFlip();
        Mesh::VertexPointer vp2dup = pos.V();
        Mesh::VertexPointer vp1 = vmap_inv[vp1dup];
        Mesh::VertexPointer vp2 = vmap_inv[vp2dup];

        // Create the two faces
        Mesh::FacePointer fp1 = &*fi++;
        Mesh::FacePointer fp2 = &*fi++;

        fp1->SetScaffold();
        fp2->SetScaffold();

        fp1->V(0) = vp2dup;
        fp1->WT(0) = pos.F()->WT(pos.VInd());

        fp1->V(1) = vp1dup;
        fp1->WT(1) = pos.F()->WT((pos.VInd()+2)%3);

        fp1->V(2) = vp1;
        fp1->WT(2) = pos.F()->WT((pos.VInd()+2)%3);

        fp2->V(0) = vp2dup;
        fp2->WT(0) = pos.F()->WT(pos.VInd());

        fp2->V(1) = vp1;
        fp2->WT(2) = pos.F()->WT((pos.VInd()+2)%3);

        fp2->V(2) = vp2;
        fp2->WT(2) = pos.F()->WT(pos.VInd());

    }

    ensure_condition(fi == m.face.end());

    /*

    for (auto pos : track) {
        Mesh::VertexPointer vp1 = pos.VFlip();
        Mesh::VertexPointer vp2 = pos.V();

        if (vmap.count(vp1) == 0) {
            vi->ImportData(*vp1);
            vmap[vp1] = &*vi;
            vi++;
        }
        if (vmap.count(vp2) == 0) {
            vi->ImportData(*vp2);
            vmap[vp2] = &*vi;
            vi++;
        }

        Mesh::VertexPointer vp1dup = vmap[pos.VFlip()];
        Mesh::VertexPointer vp2dup = vmap[pos.V()];

        // Create the two faces
        Mesh::FacePointer fp1 = &*fi++;
        Mesh::FacePointer fp2 = &*fi++;

        fp1->SetScaffold();
        fp2->SetScaffold();

        fp1->V(0) = vp2;
        fp1->WT(0) = pos.F()->WT(pos.VInd());

        fp1->V(1) = vp1;
        fp1->WT(1) = pos.F()->WT((pos.VInd()+2)%3);

        fp1->V(2) = vp1dup;
        fp1->WT(2) = pos.F()->WT((pos.VInd()+2)%3);

        fp2->V(0) = vp2;
        fp2->WT(0) = pos.F()->WT(pos.VInd());

        fp2->V(1) = vp1dup;
        fp2->WT(2) = pos.F()->WT((pos.VInd()+2)%3);

        fp2->V(2) = vp2dup;
        fp2->WT(2) = pos.F()->WT(pos.VInd());
    }
    */

    LOG_INFO << "... done.";

    return std::make_pair(vpu, fpu);
}

void AugmentMeshWithSupportFacesAlongSeams(Mesh& m)
{
    LOG_INFO << "AugmentMeshWithBlablabla";
    MarkSeamsAsFaux(m);
    tri::UpdateTopology<Mesh>::FaceFace(m);
    tri::UpdateFlags<Mesh>::FaceClearV(m);

    tri::UpdateNormal<Mesh>::PerVertex(m);

    /* Da rifare perche' cosi si invalidano gli iteratori
     * prima estraggo tutte le tracce
     * poi alloco quanto serve
     * alla fine itero su ogni traccia (facendo l'update dei puntatori)
     * */

    std::vector<std::vector<PosF>> tracks;
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            if (!f.IsV() && f.IsF(i)) {
                std::vector<PosF> fauxTrack = ExtractFauxTrack(m, PosF(&f, i));
                tracks.push_back(fauxTrack);
                for (auto& p : fauxTrack)
                    p.F()->SetV();
            }
        }
    }

    // I need to keep the pointers updated at every iteration...
    for (unsigned i = 0; i < tracks.size(); ++i) {
        std::pair<VPU,FPU> updaterPair = PadFauxTrack(m, tracks[i]);

        tri::UpdateTopology<Mesh>::FaceFace(m);
        tri::UpdateTopology<Mesh>::VertexFace(m);

        for (unsigned j = i + 1; j < tracks.size(); ++j) {
            for (auto& pos : tracks[j]) {
                updaterPair.first.Update(pos.v);
                updaterPair.second.Update(pos.f);
            }
        }
    }

    LOG_INFO << "Finished padding";

    tri::UpdateTopology<Mesh>::FaceFace(m);
    tri::UpdateTopology<Mesh>::VertexFace(m);

    for (auto& v : m.vert) {
        bool nomesh = true;
        for (face::VFIterator<Mesh::FaceType> vfi(&v); !vfi.End(); ++vfi) {
            if (vfi.F()->IsMesh()) {
                nomesh = false;
                break;
            }
        }

        if (nomesh) {
            v.P() += v.N().normalized() * 0.1;
        }
    }

    tri::io::Exporter<Mesh>::Save(m, "support.obj", tri::io::Mask::IOM_NONE);

    LOG_INFO << "Finished saving";
}





