/* This file contains classes that implement different atlas clustering
 * algorithms. All classes conform to the same interface that models a generic
 * clustering procedure as a sequence of 'moves' performed in a greedy fashion.
 * Each move, if feasible, can result in the merge of two charts into a single
 * one. The implementation of a clustering algorithm must define the following
 * components and functions.
 *
 *  - A component that computes the cost of a move, which is used to drive the
 *    greedy procedure. The lower the cost, the better the move. If a move has
 *    infinite cost, it is equivalent to be deemed unfeasbile and it cannot be
 *    performed.
 *
 *  - A component that performs the move (assuming it is feasible). This
 *    operation can have side effects on the original mesh, but must be
 *    revertible until it is declared as accepted.
 *
 *  - A component to check if the result of a merge operation is acceptable.
 *    This is a simple function that returns a boolean value.
 *
 *  - A component to finalize the move. This component either accepts the move,
 *    or rejects it. If a move is rejected, _all_ the side effects that may have
 *    occurred during the move evaluation will be reverted. Note that even if a
 *    move is rejected, the algorithm may choose to reinsert it in the queue and
 *    attempt to perform it differently at a later time (example: the algorithm
 *    detects that merging the charts along the full shared boundary does not
 *    yield satisfactory results, but can retry by mergin only a portion of the
 *    shared boundary).
 * */

#include "atlas_clustering.h"

#include "mesh_graph.h"
#include "timer.h"
#include "texture_optimization.h"
#include "matching.h"
#include "arap.h"

#include <vcg/space/index/grid_util2d.h>
#include <vcg/space/segment2.h>
#include <vcg/space/intersection2.h>
#include <vcg/complex/algorithms/attribute_seam.h>

#include <wrap/io_trimesh/export.h>


//#define FAST_CHAIN_UPDATE
//#define DEBUG_SAVES
//#define RASTERIZATION_BASED_OVERLAP_CHECK


static constexpr double OFFSET_TOLERANCE = 5.0;

static int feasibility_failed = 0;
static int post_failed = 0;
static int total_accepted = 0;
static int total_rejected = 0;

#ifdef DEBUG_SAVES

#include <vcg/complex/append.h>

static Mesh accept;
static Mesh init;

static Mesh alignmesh;
static Mesh vertmesh;

void DumpDebugMeshes(int i)
{
    std::string name = std::to_string(i) + "_accept.ply";
    tri::io::Exporter<Mesh>::Save(accept, name.c_str(), tri::io::Mask::IOM_WEDGTEXCOORD | tri::io::Mask::IOM_VERTFLAGS);

    std::string name2 = std::to_string(i) + "_init.ply";
    tri::io::Exporter<Mesh>::Save(init, name2.c_str(), tri::io::Mask::IOM_WEDGTEXCOORD | tri::io::Mask::IOM_VERTFLAGS);

    std::string name3 = std::to_string(i) + "_align.ply";
    tri::io::Exporter<Mesh>::Save(alignmesh, name3.c_str(), tri::io::Mask::IOM_WEDGTEXCOORD | tri::io::Mask::IOM_VERTFLAGS);

    std::string name4 = std::to_string(i) + "_vert.ply";
    tri::io::Exporter<Mesh>::Save(vertmesh, name4.c_str(), tri::io::Mask::IOM_VERTCOORD | tri::io::Mask::IOM_VERTCOLOR);
}
#endif


/* Remove from a boundary chain pairs of edges that cannot be merged togerher.
 * This happens when a vertex of the boundary chain is split on one side of
 * the chain but not on the other (which means that one of the two charts is
 * cut while the other is not). Note that if the same vertex is split in both
 * chains then the edges can be stitched without issues */
int RegularizeChain(std::vector<PosF>& chain1, std::vector<PosF>& chain2)
{
    ensure_condition(chain1.size() > 0);
    ensure_condition(chain1.size() == chain2.size());

    std::vector<PosF> newChain1 = { chain1[0] };
    std::vector<PosF> newChain2 = { chain2[0] };
    vcg::Point2d lastUv1 = chain1[0].F()->WT(chain1[0].VInd()).P();
    vcg::Point2d lastUv2 = chain2[0].F()->WT(chain2[0].VInd()).P();
    for (unsigned i = 1; i < chain1.size(); ++i) {
        PosF p1 = chain1[i];
        PosF p2 = chain2[i];
        p1.VFlip();
        p2.VFlip();
        vcg::Point2d uv1 = p1.F()->WT(p1.VInd()).P();
        vcg::Point2d uv2 = p2.F()->WT(p2.VInd()).P();
        bool jump1 = uv1 != lastUv1;
        bool jump2 = uv2 != lastUv1;
        if (jump1 == jump2) {
            newChain1.push_back(chain1[i]);
            newChain2.push_back(chain2[i]);
        }
        lastUv1 = chain1[i].F()->WT(chain1[i].VInd()).P();
        lastUv2 = chain2[i].F()->WT(chain2[i].VInd()).P();
    }
    int szdiff = chain1.size() - newChain1.size();
    chain1 = newChain1;
    chain2 = newChain2;

    return szdiff;
}

void ConvertChainToPoints(const BoundaryChain& chain, std::vector<Point2d>& boundary_a, std::vector<Point2d>& boundary_b)
{
    ensure_condition(chain.pos_a.size() == chain.pos_b.size());
    boundary_a.clear();
    boundary_b.clear();
    for (unsigned i = 0; i < chain.pos_a.size(); ++i) {
        const PosF& pos_a = chain.pos_a[i];
        const PosF& pos_b = chain.pos_b[i];
        boundary_a.push_back(pos_a.F()->cWT(pos_a.VInd()).P());
        boundary_b.push_back(pos_b.F()->cWT(pos_b.VInd()).P());
    }

    // Make sure at least two points are included
    if (boundary_a.size() == 1) {
        PosF pos_a = chain.pos_a[0];
        pos_a.FlipV();
        boundary_a.push_back(pos_a.F()->cWT(pos_a.VInd()).P());

        PosF pos_b = chain.pos_b[0];
        pos_b.FlipV();
        boundary_b.push_back(pos_b.F()->cWT(pos_b.VInd()).P());
    }
}

BoundaryChain ExtractBoundaryChain(Mesh& m, ChartHandle a, ChartHandle b)
{
    BoundaryChain chain;
    chain.a = a;
    chain.b = b;

    std::set<Mesh::FacePointer> visited;
    auto ccid = GetConnectedComponentIDAttribute(m);

    auto IsSeam = [&] (PosF pos) -> bool {
        auto id = ccid[pos.F()];
        auto id2 = ccid[pos.FFlip()];
        //LOG_DEBUG << "IDs are " << id << " and " << id2;
        //return ccid[pos.FFlip()] == id2;
        return (id == a->id && id2 == b->id) || (id == b->id && id2 == a->id);
    };

    PosF start;
    for (auto fptr : a->fpVec) {
        bool found = false;
        for (int i = 0; i < 3; ++i) {
            PosF p(fptr, i);
            if (IsSeam(p)) {
                start = p;
                found = true;
                break;
            }
        }
        if (found)
            break;
    }

    /* I use the fact that seams are marked as faux in the main mesh. The idea
     * of the following snippet is to start from a seam (faux) pos and walk the
     * boundary shared between charts a and b, recording pairs of edges that are
     * marked as seams (faux) in the original mesh and separate the two regions
     * _a_ and _b_ passed as arguments. What is obtained is a 'boundary chain',
     * ie a pair of sequences of topologically sorted PosF objects that
     * represents the seam.
     *
     * Note that the chain is not necessarily continuous (cuts, holes etc...) */
    PosF pos = start;
    ensure_condition(pos.IsFaux());
    ensure_condition(IsSeam(pos));


    while (true) {
        ensure_condition(pos.IsFaux());

        // If it is a seam pos, update the chain
        if (IsSeam(pos)) {
            PosF pp = pos;
            pp.FlipF();
            chain.pos_a.push_back(pos);
            chain.pos_b.push_back(pp);
        }

        PosF sp = pos;
        pos.NextFaux();

        ensure_condition(pos.IsFaux()); // ensure we remain on a faux edge
        ensure_condition(ccid[pos.F()] == ccid[start.F()]); // ensure we remain inside the correct region
        ensure_condition(pos != sp); // The pos didn't advance, something is wrong

        if (pos == start)
            break;
    }

    return chain;
}

/* Stitch a BoundaryChain by merging the uv coordinates along edges. This
 * function iterates over each pair of PosF objects, computes a new uv point as
 * the average of the two uv points on each side of the chain, then iterates
 * over the faces incident to the vertex pointed by the PosF, assigning the
 * new uv coordinate to the vertices that have the same initial uv in each side.
 * This allows to gracefully handle cuts and non-manifoldness in uv space. */
StitchOffsetVec StitchChain(Mesh& m, const BoundaryChain& chain)
{
    auto ccid = GetConnectedComponentIDAttribute(m);

    StitchOffsetVec offsetVec;

    /* First pass */
    for (unsigned i = 0; i < chain.pos_a.size(); ++i) {
        PosF pa = chain.pos_a[i];
        PosF pb = chain.pos_b[i];
        vcg::Point2d uva = pa.F()->WT(pa.VInd()).P();
        vcg::Point2d uvb = pb.F()->WT(pb.VInd()).P();
        vcg::Point2d uv = 0.5 * (uva + uvb);

        double offset = OFFSET_TOLERANCE * std::max((uva - uv).Norm(), (uvb - uv).Norm());
        offsetVec.push_back(std::make_pair(uv, offset));

        PosF pos = pa;
        do {
            if (pos.F()->WT(pos.VInd()).P() == uva && (ccid[pos.F()] == chain.a->id)) {
                pos.F()->WT(pos.VInd()).P() = uv;
            } else if (pos.F()->WT(pos.VInd()).P() == uvb && (ccid[pos.F()] == chain.b->id)) {
                pos.F()->WT(pos.VInd()).P() = uv;
            }
            pos.NextE();
        } while (pos != pa);
    }


    /* This is the ugly bit. We need to ensure that each edge has been stitched
     * along both vertices, so flip each pos and check if the coordinate is the
     * same of the previous pos (in topological order). If there is a jump on
     * both sides, then we must also stitch using the flipped pos */
    for (unsigned i = 0; i < chain.pos_a.size(); ++i) {
        PosF ppa = chain.pos_a[i];
        PosF ppb = chain.pos_b[i];
        ppa.FlipV();
        ppb.FlipV();
        unsigned j = (i > 0) ? i-1 : chain.pos_a.size() - 1;

        PosF pa = chain.pos_a[j];
        PosF pb = chain.pos_b[j];
        bool jmpa = ppa.F()->WT(ppa.VInd()).P() != pa.F()->WT(pa.VInd()).P();
        bool jmpb = ppb.F()->WT(ppb.VInd()).P() != pb.F()->WT(pb.VInd()).P();

        if (jmpa && jmpb) {
            vcg::Point2d uva = ppa.F()->WT(ppa.VInd()).P();
            vcg::Point2d uvb = ppb.F()->WT(ppb.VInd()).P();
            vcg::Point2d uv = 0.5 * (uva + uvb);

            double offset = OFFSET_TOLERANCE * std::max((uva - uv).Norm(), (uvb - uv).Norm());
            offsetVec.push_back(std::make_pair(uv, offset));

            PosF pos = ppa;
            do {
                if (pos.F()->WT(pos.VInd()).P() == uva && (ccid[pos.F()] == chain.a->id)) {
                    pos.F()->WT(pos.VInd()).P() = uv;
                } else if (pos.F()->WT(pos.VInd()).P() == uvb && (ccid[pos.F()] == chain.b->id)) {
                    pos.F()->WT(pos.VInd()).P() = uv;
                }
                pos.NextE();
            } while (pos != ppa);
        }

        // Remove the faux flag from the edge pair, since we are closing the seam
        if (jmpa == jmpb) {
            pa.F()->ClearF(pa.E());
            pb.F()->ClearF(pb.E());
        }
    }

    return offsetVec;
}


/* Atlas clustering implementation
 * =============================== */

AtlasClustering::AtlasClustering(std::shared_ptr<MeshGraph> gptr)
    : g{gptr},
      queue{},
      moves{},
      chains{},
      total_arap_t{0},
      active_arap_t{0},
      init_t{0},
      check_t{0},
      merge_t{0},
      opt_t{0},
      post_t{0},
      shell_0_t{0},
      shell_1_t{0},
      shell_2_t{0},
      shell_3_t{0}

{
    ResetMoveState();

    for (auto & elem : g->charts)
        elem.second->ParameterizationChanged();

    tri::UpdateTopology<Mesh>::FaceFace(g->mesh);
    for (auto elem : g->charts) {
        ChartHandle c1 = elem.second;
        for (ChartHandle c2 : c1->adj) {
            ClusteringMove mv(c1, c2);
            if (moves.count(mv) == 0) {
                WeightedClusteringMove wm;
                bool res = ComputeWeightedMove(mv, wm, chains[mv]);
                ensure_condition(res);
                AddMove(wm);
            }
        }
    }
}

std::pair<AtlasClustering::MergeStatus,std::queue<ChartHandle>> AtlasClustering::MergeCheck(const std::vector<ChartHandle> charts)
{
    auto first = charts.begin();
    auto last = charts.end();
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
        auto p = std::make_pair(MergeStatus::Allowed, std::queue<ChartHandle>{});
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
        return std::make_pair(MergeStatus::ErrorDisconnected, std::queue<ChartHandle>{});
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
        return std::make_pair(MergeStatus::Allowed, mergeQueue);
    } else {
        return std::make_pair(MergeStatus::ErrorUnfeasible, std::queue<ChartHandle>{});
    }
}

bool AtlasClustering::ComputeWeightedMove(const ClusteringMove& mv, WeightedClusteringMove& wm, BoundaryChain& chain)
{
    chain = ExtractBoundaryChain(g->mesh, mv.a, mv.b);
    double cost = ComputeCost(mv, chain);
    wm.first = mv;
    wm.second = cost;
    return (chain.pos_a.size() > 0);
}

double AtlasClustering::ComputeCost(const ClusteringMove& mv, const BoundaryChain& chain)
{
    ensure_condition(mv.a == chain.a && mv.b == chain.b);
    double sharedUV_a = 0;
    double sharedUV_b = 0;
    for (unsigned i = 0; i < chain.pos_a.size(); ++i) {
        sharedUV_a += EdgeLengthUV(*(chain.pos_a[i].F()), chain.pos_a[i].E());
        sharedUV_b += EdgeLengthUV(*(chain.pos_b[i].F()), chain.pos_b[i].E());
    }
    double maxbndf = std::max(sharedUV_a / mv.a->BorderUV(), sharedUV_b / mv.b->BorderUV());
    double penalty = std::max(sharedUV_a / sharedUV_b, sharedUV_b / sharedUV_a);
    double cost = std::max(1 - maxbndf, 0.0) * penalty;
    return cost;
}

bool AtlasClustering::AddMove(const WeightedClusteringMove& wm)
{
    const ClusteringMove& mv = wm.first;
    ensure_condition(mv.a->adj.count(mv.b) > 0 && mv.b->adj.count(mv.a) > 0);
    if (moves.find(mv) == moves.end()) {
        moves.insert(wm);
        queue.push(wm);
        return true;
    }
    else {
        return false;
    }
}

void AtlasClustering::DeleteMove(const ClusteringMove& move)
{
    moves.erase(move);
    chains.erase(move);
}

void AtlasClustering::InvalidateMove(const ClusteringMove& m)
{
    moves[m] = Infinity();
    queue.push(std::make_pair(m, Infinity()));
}

/*
 * Move steps implementation
 * =========================
 * */

void AtlasClustering::ResetMoveState()
{
    state = Uninitialized;
    currentMove = ClusteringMove::NullMove();
    merge_a = nullptr;
    merge_b = nullptr;
    texIdA = -1;
    texIdB = -1;
    offsetVec.clear();
    savedWtA.clear();
    savedWtB.clear();
    move_arap_t = 0;
}

void AtlasClustering::InitMove(ClusteringMove move)
{
    Timer t;

    LOG_DEBUG << "    Running  step InitMove()";
    ensure_condition(state == Uninitialized);
    currentMove = move;
    merge_a = move.a;
    merge_b = move.b;
    state = Initialized;
    LOG_DEBUG << "Initialized move " << merge_a->id << " <- " << merge_b->id;

    init_t += t.TimeElapsed();
}

bool AtlasClustering::CheckFeasibility()
{
    Timer t;

    LOG_DEBUG << "    Running step CheckFeasibility()";
    ensure_condition(state == Initialized);
    std::vector<ChartHandle> merge = { merge_a, merge_b };
    auto check = MergeCheck(merge);
    if (check.first == MergeStatus::Allowed) {
        state = Feasible;
    } else {
        LOG_DEBUG << "Unfeasible move";
        state = Unfeasible;
    }

    // Additionally, also check the alignment error
    ensure_condition(chains.count(currentMove) > 0);
    const BoundaryChain& chain = chains[currentMove];

    std::vector<vcg::Point2d> bpa;
    std::vector<vcg::Point2d> bpb;

    ConvertChainToPoints(chain, bpa, bpb);

    MatchingInfo mi = ComputeMatchingRigidMatrix(bpa, bpb);

    double error = 0;
    for (unsigned i = 0; i < bpa.size(); ++i) {
        error += (bpa[i] - mi.Apply(bpb[i])).Norm();
    }

    error = error / (double) bpa.size();

    const double ERROR_THRESHOLD = 20;

    if (error > ERROR_THRESHOLD) {
        LOG_DEBUG << "Forcing infeasibility because matching error is too high (" << error << ")";
        state = Unfeasible;
    }

    bool passed = (state == Feasible);
    if (!passed)
        feasibility_failed++;

    check_t += t.TimeElapsed();

    return passed;
}


void AtlasClustering::Merge()
{
    Timer t;

    LOG_DEBUG << "    Running step Merge()";
    ensure_condition(state == Feasible);
    /* In order to be able to revert the changes, after merging the charts
     * I need to store the reference to the 'new' charts, the number of faces
     * that have been added because of the merge, and the texture coordinates of
     * the faces of both charts. */
    savedWtA.clear();
    savedWtA.reserve(6 * merge_a->FN());
    for (auto fptr : merge_a->fpVec) {
        for (int i = 0; i < fptr->VN(); ++i) {
            savedWtA.push_back(fptr->WT(i).U());
            savedWtA.push_back(fptr->WT(i).V());
        }
    }
    texIdA = merge_a->fpVec[0]->WT(0).N();

    savedWtB.clear();
    savedWtB.reserve(6 * merge_b->FN());
    for (auto fptr : merge_b->fpVec) {
        for (int i = 0; i < fptr->VN(); ++i) {
            savedWtB.push_back(fptr->WT(i).U());
            savedWtB.push_back(fptr->WT(i).V());
        }
    }
    texIdB = merge_b->fpVec[0]->WT(0).N();

    /* Perform the merge of b into a without changing the adjacencies. Note that
     * this leaves the mesh graph in an inconsistent state. Maybe it would be
     * better to perform the merge and then revert it by using the reference to
     * b, which has all the correct references to the other atlas charts
     * still... */

    ensure_condition(chains.count(currentMove) > 0);
    const BoundaryChain& chain = chains[currentMove];

    std::vector<vcg::Point2d> bpa;
    std::vector<vcg::Point2d> bpb;

    ConvertChainToPoints(chain, bpa, bpb);

    MatchingInfo mi = ComputeMatchingRigidMatrix(bpa, bpb);

    LOG_DEBUG << "Matching info";
    LOG_DEBUG << mi.t.X() << " " << mi.t.Y();
    LOG_DEBUG << mi.matCoeff[0] << " " << mi.matCoeff[1] << " " << mi.matCoeff[2] << " " << mi.matCoeff[3];

    double error = 0;
    for (unsigned i = 0; i < bpa.size(); ++i) {
        error += (bpa[i] - mi.Apply(bpb[i])).Norm();
    }

    LOG_DEBUG << "Alignment error is " << error;
    merge_a->error += error;

    // Align merge_b to merge_a, and change the texture id to match
    for (auto fptr : merge_b->fpVec) {
        for (int i = 0; i < fptr->VN(); ++i) {
            fptr->WT(i).P() = mi.Apply(fptr->WT(i).P());
            fptr->WT(i).N() = texIdA;
        }
    }

#ifdef DEBUG_SAVES
    {
        //static int i = 0;
        alignmesh.Clear();
        alignmesh.ClearAttributes();
        Mesh m2;
        MeshFromFacePointers(merge_a->fpVec, alignmesh);
        MeshFromFacePointers(merge_b->fpVec, m2);
        tri::Append<Mesh,Mesh>::Mesh(alignmesh, m2);

        for (auto& f : alignmesh.face) {
            for (int i = 0; i < 3; ++i) {
                f.P(i) = Point3d(f.WT(i).U(), f.WT(i).V(), 0);
            }
        }

        Mesh ma;
        CopyToMesh(*merge_a, ma);
        Mesh mb;
        CopyToMesh(*merge_b, mb);

        tri::io::Exporter<Mesh>::Save(ma, "ma.ply", tri::io::Mask::IOM_WEDGTEXCOORD | tri::io::Mask::IOM_VERTFLAGS);
        tri::io::Exporter<Mesh>::Save(mb, "mb.ply", tri::io::Mask::IOM_WEDGTEXCOORD | tri::io::Mask::IOM_VERTFLAGS);


        vertmesh.Clear();
        for (auto p : bpa) {
            auto vi = tri::Allocator<Mesh>::AddVertex(vertmesh, vcg::Point3d(p.X(), p.Y(), 0));
            vi->C() = vcg::Color4b::Red;
        }


        for (auto pp : bpb) {
            vcg::Point2d p = mi.Apply(pp);
            auto vi = tri::Allocator<Mesh>::AddVertex(vertmesh, vcg::Point3d(p.X(), p.Y(), 0));
            vi->C() = vcg::Color4b::Blue;
        }
    }
#endif

    // Stitch the charts along the active boundary chain
    offsetVec = StitchChain(g->mesh, chain);

    merge_a->fpVec.reserve(merge_a->fpVec.size() + merge_b->fpVec.size());
    for (auto fptr : merge_b->fpVec)
        merge_a->AddFace(fptr);

    state = Merged;

    merge_t += t.TimeElapsed();

    OptimizeAfterMerge();
}

void AtlasClustering::OptimizeAfterMerge()
{
    Timer t;

    LOG_DEBUG << "    Running step OptimizeAfterMerge()";

    ensure_condition(offsetVec.size() > 0);
    ensure_condition(state == MoveState::Merged);

    LOG_DEBUG << "Performing optimization after stitching...";

    LOG_DEBUG << "Building shell...";

    Mesh shell;
    //shell.Clear();
    //shell.ClearAttributes();
    BuildShell(shell, *merge_a, ParameterizationGeometry::Texture, false);
    ClearHoleFillingFaces(shell, true, true);

    shell_0_t += t.TimeElapsed();

    // Copy the texture coordinates per wedge from the clustered chart
    for (unsigned i = 0; i < merge_a->FN(); ++i) {
        auto& sf = shell.face[i];
        auto& f = *(merge_a->fpVec[i]);
        for (int j = 0; j < 3; ++j)
            sf.WT(j) = f.WT(j);
    }

    // Cut along seams (convert per wedge uvs into per vertex uvs)
    auto vExt = [](const Mesh& msrc, const MeshFace& f, int k, const Mesh& mdst, MeshVertex& v) {
        (void) msrc;
        (void) mdst;
        v.ImportData(*(f.cV(k)));
        v.T() = f.cWT(k);
    };
    auto vCmp = [](const Mesh& mdst, const MeshVertex& v1, const MeshVertex& v2) {
        (void) mdst;
        return v1.T() == v2.T();
    };

    int vn = shell.vn;

    shell_1_t += t.TimeElapsed();

    tri::AttributeSeam::SplitVertex(shell, vExt, vCmp);
    if (shell.VN() != vn) {
        tri::Allocator<Mesh>::CompactEveryVector(shell);
        tri::UpdateTopology<Mesh>::FaceFace(shell);
        tri::UpdateTopology<Mesh>::VertexFace(shell);
    }

    SyncShellWithUV(shell);

    /*
    for (auto& sf : shell.face)
        for (int i = 0; i < 3; ++i)
            sf.V(i)->T() = sf.WT(i);
            */


    shell_2_t += t.TimeElapsed();

#ifdef DEBUG_SAVES
    init.Clear();
    init.ClearAttributes();
    tri::Append<Mesh, Mesh>::MeshCopy(init, shell);
#endif


    // Restrict the shell, leave only the faces that have at least one texture coordinate
    // within the tolerance of any stitched uv point

    tri::UpdateFlags<Mesh>::VertexSetS(shell);
    int del = 0;
    int cleared = 0;
    for (auto& sf : shell.face) {
        bool inside[3] = { false, false, false };
        for (int i = 0; i < 3; ++i) {
            inside[i] = !sf.V(i)->IsS();
            if (!inside[i]) {
                for (const auto& entry : offsetVec) {
                    if ((sf.WT(i).P() - entry.first).Norm() <= entry.second) {
                        inside[i] = true;
                        sf.V(i)->ClearS();
                        cleared++;
                        break;
                    }
                }
            }
        }
        if (!inside[0] && !inside[1] && !inside[2]) {
            tri::Allocator<Mesh>::DeleteFace(shell, sf);
            del++;
        }
    }

    if (del > 0) {
        LOG_DEBUG << "Shell restriction deleted " << del << " faces";
        tri::Clean<Mesh>::RemoveUnreferencedVertex(shell);
        tri::Allocator<Mesh>::CompactEveryVector(shell);
        tri::UpdateTopology<Mesh>::FaceFace(shell);
        tri::UpdateTopology<Mesh>::VertexFace(shell);
    }

    shell_3_t += t.TimeElapsed();

    Timer arap_timer;

    LOG_DEBUG << "Building ARAP object...";
    ARAP arap(shell);
    arap.SetMaxIterations(50);

    // Fix the selected vertices of the shell;
    LOG_DEBUG << "Fixing selected vertices";
    int nfixed = arap.FixSelectedVertices();
    ensure_condition(nfixed == (shell.VN() - cleared));
    double tol = 0.02;
    while (nfixed < 2) {
        LOG_DEBUG << "Not enough selected vertices found, fixing random edge with tolerance " << tol;
        nfixed += arap.FixRandomEdgeWithinTolerance(tol);
        tol += 0.02;
    }
    ensure_condition(nfixed > 0);

    // Fix the boundary of the shell
    //tri::UpdateFlags<Mesh>::VertexBorderFromFaceAdj(shell);
    //arap.FixBoundaryVertices();

    LOG_DEBUG << "Solving...";
    arap.Solve();

    move_arap_t = arap_timer.TimeElapsed();
    total_arap_t += move_arap_t;

    SyncShellWithUV(shell);

    LOG_DEBUG << "Syncing chart...";

    ensure_condition(HasFaceIndexAttribute(shell));
    auto ia = GetFaceIndexAttribute(shell);
    for (auto& sf : shell.face) {
        auto& f = (g->mesh).face[ia[sf]];
        for (int k = 0; k < 3; ++k)
            f.WT(k).P() = sf.V(k)->T().P();
    }

#ifdef DEBUG_SAVES
    accept.Clear();
    accept.ClearAttributes();

    BuildShell(accept, *merge_a, ParameterizationGeometry::Texture, false);
    ClearHoleFillingFaces(accept, true, true);

    int vn = shell.VN();
    tri::AttributeSeam::SplitVertex(accept, vExt, vCmp);
    if (accept.VN() != vn) {
        tri::Allocator<Mesh>::CompactEveryVector(accept);
        tri::UpdateTopology<Mesh>::FaceFace(accept);
        tri::UpdateTopology<Mesh>::VertexFace(accept);
    }

    /*
    for (auto& sf : shell.face)
        for (int i = 0; i < 3; ++i)
            sf.V(i)->T() = sf.WT(i);
            */

    SyncShellWithUV(accept);

    //tri::Append<Mesh, Mesh>::MeshCopy(accept, shell);
#endif


    opt_t += t.TimeElapsed();
}

static bool SegmentBoxIntersection(const Segment2<double>& seg, const Box2d& box)
{
    Point2d isec;
    Point2d c1{box.min};
    Point2d c2{box.max[0], box.min[1]};
    Point2d c3{box.max};
    Point2d c4{box.min[0], box.max[1]};

    if (SegmentSegmentIntersection(seg, Segment2<double>{c1, c2}, isec)) return true;
    if (SegmentSegmentIntersection(seg, Segment2<double>{c2, c3}, isec)) return true;
    if (SegmentSegmentIntersection(seg, Segment2<double>{c3, c4}, isec)) return true;
    if (SegmentSegmentIntersection(seg, Segment2<double>{c4, c1}, isec)) return true;

    // if the segment does not intersect the sides, check if it is fully contained in the box
    return (box.min[0] <= std::min(seg.P0()[0], seg.P1()[0]) &&
            box.min[1] <= std::min(seg.P0()[1], seg.P1()[1]) &&
            box.max[0] >= std::max(seg.P0()[0], seg.P1()[0]) &&
            box.max[1] >= std::max(seg.P0()[1], seg.P1()[1]));
}

#include "texture_rendering.h"
static bool OverlapCheck(ChartHandle merge_a, double frac)
{
    RasterizedParameterizationStats stats = GetRasterizationStats(merge_a, 1024, 1024);
    double fraction = stats.lostFragments / (double) stats.totalFragments;
    if (fraction > frac) {
        LOG_VERBOSE << "OverlapCheck() failed, overlap fraction = " << fraction << ")";
        return false;
    }
    return true;
}

#ifdef RASTERIZATION_BASED_OVERLAP_CHECK
bool AtlasClustering::PostCheck()
{
    Timer t;
    LOG_DEBUG << "    Running step PostCheck()";
    ensure_condition(state == Merged);
    state = Acceptable;

    const double fail_threshold = 0.005;

    RasterizedParameterizationStats stats = GetRasterizationStats(merge_a, 1024, 1024);
    double fraction = stats.lostFragments / (double) stats.totalFragments;
    if (fraction > fail_threshold) {
        LOG_VERBOSE << "PostCheck() failed, overlap fraction = " << fraction << ")";
        post_failed++;
    }

    post_t += t.TimeElapsed();

    return !(fraction > fail_threshold);
}
#else
bool AtlasClustering::PostCheck()
{
    Timer t;

    LOG_DEBUG << "    Running step PostCheck()";
    ensure_condition(state == Merged);
    state = Acceptable;

    MyMesh em1;
    MyMesh em2;

    /* retrieve a collection of non-seam edges using the fact that the ids of the
     * faces of merge_a have not been updated yet. Then check using a hashgrid
     * if there are any edges from the original merge_b face that overlap with
     * the boundary of merge_a. Note that this is not sufficient in the case
     * where we merge only a sub-portion of the boundary edges, but we can
     * resort to additionally check if the shared edges has been stitched. */

    using EdgeUV = vcg::Segment2<double>;

    /* edges[0] contains the edges of merge_a that were not stitched
     * edges[1] contains the edges of merge_b that were not stitched */
    std::vector<std::vector<EdgeUV>> edges(2);
    auto ccid = GetConnectedComponentIDAttribute(g->mesh);
    for (auto fptr : merge_a->fpVec) {
        auto idf = ccid[fptr];
        for (int i = 0; i < 3; ++i) {
            if (IsBorderUV(fptr, i)) {
                int k = (idf == merge_a->id) ? 0 : 1;
                edges[k].push_back(EdgeUV(fptr->WT(i).P(), fptr->WT((i+1)%3).P()));
                vcg::Point2d p1 = fptr->WT(i).P();
                vcg::Point2d p2 = fptr->WT((i+1)%3).P();
                if (k == 0)
                    tri::Allocator<MyMesh>::AddEdge(em1, vcg::Point3d(p1.X(), p1.Y(), 0), vcg::Point3d(p2.X(), p2.Y(), 0));
                else
                    tri::Allocator<MyMesh>::AddEdge(em2, vcg::Point3d(p1.X(), p1.Y(), 0), vcg::Point3d(p2.X(), p2.Y(), 0));
            }
        }
    }

    /*
    bool check = false;
    for (auto fptr : merge_a->fpVec) {
        if (tri::Index(g->mesh, fptr) == 144098)
            check = true;
    }

    if (check) {
        static int i = 0;
        std::string n1 = "em1_" + std::to_string(i) + ".ply";
        std::string n2 = "em2_" + std::to_string(i) + ".ply";
        tri::io::Exporter<MyMesh>::Save(em1, n1.c_str(), tri::io::Mask::IOM_VERTCOORD | tri::io::Mask::IOM_EDGEINDEX);
        tri::io::Exporter<MyMesh>::Save(em2, n2.c_str(), tri::io::Mask::IOM_VERTCOORD | tri::io::Mask::IOM_EDGEINDEX);

        auto vExt = [](const Mesh& msrc, const MeshFace& f, int k, const Mesh& mdst, MeshVertex& v) {
            (void) msrc;
            (void) mdst;
            v.ImportData(*(f.cV(k)));
            v.T() = f.cWT(k);
        };
        auto vCmp = [](const Mesh& mdst, const MeshVertex& v1, const MeshVertex& v2) {
            (void) mdst;
            return v1.T() == v2.T();
        };

        Mesh shell;
        BuildShell(shell, *merge_a, ParameterizationGeometry::Texture, false);
        ClearHoleFillingFaces(shell, true, true);

        int vn = shell.VN();
        tri::AttributeSeam::SplitVertex(shell, vExt, vCmp);
        if (shell.VN() != vn) {
            tri::Allocator<Mesh>::CompactEveryVector(shell);
            tri::UpdateTopology<Mesh>::FaceFace(shell);
            tri::UpdateTopology<Mesh>::VertexFace(shell);
        }

        SyncShellWithUV(shell);

        std::string sn = "shell_" + std::to_string(i) + ".ply";
        tri::io::Exporter<Mesh>::Save(shell, sn.c_str(), tri::io::Mask::IOM_VERTCOORD);

        i++;
    }
    */

    using SegmentID = std::pair<int, int>;

    std::unordered_map<vcg::Point2i, std::vector<SegmentID>, Point2iHasher> grid;
    unsigned elems = edges[0].size() + edges[1].size();

    // init grid helper
    vcg::BasicGrid2D<double> gh;
    gh.bbox = merge_a->UVBox();
    BestDim2D<double>(elems, gh.bbox.Dim(), gh.siz);
    gh.ComputeDimAndVoxel();

    for (int i = 0; i < 2; ++i) {
        for (unsigned j = 0; j < edges[i].size(); ++j) {
            const EdgeUV& e = edges[i][j];

            Box2d edgeBox;
            edgeBox.Add(e.P0());
            edgeBox.Add(e.P1());
            Box2i gridCover;
            gh.BoxToIBox(edgeBox, gridCover);
            for (int h = gridCover.min[0]; h <= gridCover.max[0]; h++) {
                for (int k = gridCover.min[1]; k <= gridCover.max[1]; k++) {
                    Box2d cell;
                    Point2i voxel{h, k};
                    gh.IPiToBox(voxel, cell);
                    if (SegmentBoxIntersection(e, cell)) {
                        grid[voxel].push_back(std::make_pair(i, j));
                    }
                }
            }
        }
    }

    for (auto& entry : grid) {
        for (auto i1 : entry.second) {
            for (auto i2 : entry.second) {
                /* skip test if the elements are the same */
                if (i1 == i2)
                    continue;

                /* skip the test also if the edges belong to the same chart,
                 * this way we avoid rejecting due to overlaps that were present
                 * BEFORE the merge took place. Note however that this is only
                 * valid if the post-merge optimization leaves the boundary
                 * unchanged (which as of now is not guaranteed...) */
                if (i1.first == i2.first)
                    continue;

                Segment2<double> s1 = edges[i1.first][i1.second];
                Segment2<double> s2 = edges[i2.first][i2.second];
                Point2d intersectionPoint;
                if (SegmentSegmentIntersection(s1, s2, intersectionPoint)
                        && intersectionPoint != s1.P0() && intersectionPoint != s1.P1()
                        && intersectionPoint != s2.P0() && intersectionPoint != s2.P1()) {
                    LOG_DEBUG << "( " <<  s1.P0()[0] << " , " << s1.P0()[1] << " ) (" << s1.P1()[0] << " , " << s1.P1()[1] << " )";
                    LOG_DEBUG << "( " <<  s2.P0()[0] << " , " << s2.P0()[1] << " ) (" << s2.P1()[0] << " , " << s2.P1()[1] << " )";
                    LOG_DEBUG << intersectionPoint[0] << " , " << intersectionPoint[1];
                    post_failed++;
                    post_t += t.TimeElapsed();
                    return false;
                }
            }
        }
    }

    post_t += t.TimeElapsed();
    return true;
}
#endif

/* Notes on what we need to update when accepting a move.
 * Essentially, when we accept a move we change the set of valid moves; however,
 * moves are deeply linked with their boundary chain. There are three different
 * scenarios in which a chart can find itself:
 *   1. A chart c was adjacent to a, but not to b
 *   2. A chart c was adjacent to b, but not to a,
 *   3. A chart c was adjacent to both a and b
 *
 * (1) The BoundaryChain between a and c does not change, but the cost must
 *     be updated
 * (2) The BoundaryChain between a and c is the same as the boundary chain
 *     between b and c (unless it was shortened, TODO FIXME). In this case
 *     the chain must be copied, and the cost must be updated/recomputed
 * (3) In this case the chain between a and c is not maximal, we simply
 *     readd the move entirely
 * */
void AtlasClustering::AcceptMove()
{
    LOG_DEBUG << "    Running step AcceptMove()";
    using ChartSet = std::unordered_set<ChartHandle, FaceGroup::Hasher>;

    ensure_condition(state == Acceptable || state == Unacceptable);

    auto ccid = GetConnectedComponentIDAttribute(g->mesh);

    // Update ID and add faces
    for (auto fptr : merge_b->fpVec)
        ccid[fptr] = merge_a->id;

#ifdef FAST_CHAIN_UPDATE
    // Trying to be clever here....
    ChartSet shared;
    ChartSet only_a = merge_a->adj;
    ChartSet only_b = merge_b->adj;

    for (auto ch : merge_a->adj) {
        if (only_b.erase(ch) > 0) {
            int n = only_a.erase(ch);
            ensure_condition(n > 0);
            shared.insert(ch);
        }
    }
    ensure_condition(only_a.size() + only_b.size() + 2 * shared.size() == merge_a->adj.size() + merge_b->adj.size());
#endif

    // Update adjacencies in the atlas graph
    merge_a->adj.erase(merge_b);
    for (auto c : merge_b->adj) {
        if (c != merge_a) { // merge_a is now (if it wasn't already) adjacent to c
            c->adj.erase(merge_b);
            c->adj.insert(merge_a);
            merge_a->adj.insert(c);
        }
    }

#ifdef FAST_CHAIN_UPDATE
    for (auto ch : merge_a->adj) {
        ClusteringMove mv(merge_a, ch);
        if (shared.count(ch)) {
            // Recompute the move entirely
            DeleteMove(mv);
            WeightedClusteringMove wm;
            bool res = ComputeWeightedMove(mv, wm, chains[mv]);
            ensure_condition(res);
            AddMove(wm);
        } else if (only_a.count(ch)) {
            // Boundary unchanged, only update the cost of the move
            double cost = ComputeCost(mv, chains[mv]);
            DeleteMove(mv);
            AddMove(std::make_pair(mv, cost));
        } else if (only_b.count(ch)) {
            // Re-use the old chain instead of recomputing it
            ClusteringMove oldMove(merge_b, ch);
            ensure_condition(chains.count(oldMove) > 0);
            chains[mv] = chains[oldMove];
            //BoundaryChain& oldChain = chains[oldMove];

            // Update chart reference in the old chain
            if (chains[mv].a == merge_b) {
                chains[mv].a = merge_a;
            } else {
                ensure_condition(chains[mv].b == merge_b);
                chains[mv].b = merge_a;
            }

            // Swap references in the old chain to conform to the order of mv
            if (chains[mv].a == mv.b) {
                ensure_condition(chains[mv].b == mv.a);
                std::swap(chains[mv].a, chains[mv].b);
                std::swap(chains[mv].pos_a, chains[mv].pos_b);
            }

            // Sanity checks
            ensure_condition(chains[mv].a == mv.a && chains[mv].b == mv.b);
            ensure_condition(chains[mv].a->id == ccid[chains[mv].pos_a[0].F()]);
            ensure_condition(chains[mv].b->id == ccid[chains[mv].pos_b[0].F()]);

            chains.erase(oldMove);

            /*

            if (oldChain.a == merge_b)
                oldChain.a = merge_a;
            else
                oldChain.a = merge_b;

            // make sure the chain is consistent with the chart refs of the move
            if (merge_a == oldChain.b) {
                std::swap(oldChain.a, oldChain.b);
                std::swap(oldChain.pos_a, oldChain.pos_b);
            }
            */

            // Update cost of the move
            double cost = ComputeCost(mv, chains[mv]);
            DeleteMove(mv);
            AddMove(std::make_pair(mv, cost));
        } else {
            // Should never happen if shared and only_* are computed correctly
            ensure_condition(0 && "Error with adjacent charts partitioning");
        }

        DeleteMove(ClusteringMove(merge_b, ch));
    }

#endif

    for (auto ch : merge_a->adj) {
        ClusteringMove mv(merge_a, ch);

        // Recompute the move entirely
        DeleteMove(mv);
        WeightedClusteringMove wm;
        bool res = ComputeWeightedMove(mv, wm, chains[mv]);
        ensure_condition(res);
        AddMove(wm);

        DeleteMove(ClusteringMove(merge_b, ch));
    }


    merge_a->numMerges += merge_b->numMerges + 1;

    g->charts.erase(merge_b->id);

    total_accepted++;
    active_arap_t += move_arap_t;
    LOG_DEBUG << "Accepted move " << merge_a->id << " <- " << merge_b->id << " (" << total_accepted << " total)";

#ifdef DEBUG_SAVES
    DumpDebugMeshes(total_accepted);
#endif

    moves.erase(currentMove);

    ResetMoveState();
}

void AtlasClustering::RejectMove()
{
    LOG_DEBUG << "    Running step RejectMove()";
    ensure_condition(state == Acceptable || state == Unacceptable || state == Unfeasible);

    /* To revert the move:
     *   Assuming n faces have been added from b to a:
     *     remove the last n faces from a->fpVec
     *     restore the uv coordinates (if they have changed) for both charts
     * Then possibly recompute a less aggressive move, and put it back in the
     * queue.
     * Finally, update the state of the object */
    if (state == Acceptable || state == Unacceptable) {

#ifdef DEBUG_SAVES
        if (savedWtA.size() == 6) {
            LOG_DEBUG << savedWtA[0] << " , " << savedWtA[1];
            LOG_DEBUG << savedWtA[2] << " , " << savedWtA[3];
            LOG_DEBUG << savedWtA[4] << " , " << savedWtA[5];

        }

        if (savedWtB.size() == 6) {
            LOG_DEBUG << savedWtB[0] << " , " << savedWtB[1];
            LOG_DEBUG << savedWtB[2] << " , " << savedWtB[3];
            LOG_DEBUG << savedWtB[4] << " , " << savedWtB[5];

        }

        if (savedWtA.size() == 6 || savedWtB.size() == 6) {
            DumpDebugMeshes(total_rejected);
            assert(0);
        }
#endif



        merge_a->fpVec.erase(merge_a->fpVec.end() - merge_b->FN(), merge_a->fpVec.end());
        double *uvptra = savedWtA.data();
        for (auto fptr : merge_a->fpVec) {
            for (int i = 0; i < fptr->VN(); ++i) {
                fptr->WT(i).U() = *uvptra++;
                fptr->WT(i).V() = *uvptra++;
            }
        }
        merge_a->ParameterizationChanged();
        double *uvptrb = savedWtB.data();
        for (auto fptr : merge_b->fpVec) {
            for (int i = 0; i < fptr->VN(); ++i) {
                fptr->WT(i).U() = *uvptrb++;
                fptr->WT(i).V() = *uvptrb++;
                fptr->WT(i).N() = texIdB;
            }
        }
        merge_b->ParameterizationChanged();

        // Restore the faux flag for all the edges that are in the chain (since
        // the merge is rejected they remain seam edges)
        BoundaryChain& chain = chains[currentMove];
        for (unsigned i = 0; i < chain.pos_a.size(); ++i) {
            PosF& pa = chain.pos_a[i];
            PosF& pb = chain.pos_b[i];
            pa.F()->SetF(pa.E());
            pb.F()->SetF(pb.E());
        }

        // TODO: here we could manipulate the extra data linked to the move (the boundary
        // chain that we want to stitch together
        /* example:
         *  chains[currentMove] = UpdateChain(chains[currentMove]);
         *  UpdateCost(currentMove)
         * */
        // for now, just invalidate the move
        InvalidateMove(currentMove);
    } else {
        moves.erase(currentMove);
    }

    ResetMoveState();
    total_rejected++;
}


/* The clustering loop */
void AtlasClustering::Run(std::size_t targetAtlasSize, double smallRegionAreaThreshold)
{
    Timer timer;

    ensure_condition(targetAtlasSize > 0);

    LOG_INFO << "Atlas clustering: target size is " << targetAtlasSize << ", small area threshold is " << smallRegionAreaThreshold;

    double minChartArea = smallRegionAreaThreshold * g->Area3D();


    /*
    LOG_DEBUG << "State dump";
    LOG_DEBUG << "Number of nodes: " << g->charts.size() << " (node list follows)";
    for (auto& entry : g->charts) {
        LOG_DEBUG << "Chart " << entry.first << " (FN = " << entry.second->FN() << ")";
    }
    LOG_DEBUG << "Number of moves: " << moves.size() << " (edge list follows)";
    for (const auto& entry : moves) {
        LOG_DEBUG << entry.first.a->id << " " << entry.first.b->id << " | w = " << entry.second;
    }
    */

    int mergeCount;
    int moveCount;
    int numIter = 0;
    do {
        //mergeCount = CloseMacroRegions(smallRegionAreaThreshold);

        mergeCount = 0;
        moveCount = 0;

        while (InitializeNextMove()) {

            ClusteringMove move = CurrentMove();
            double minNextArea = std::min(move.a->Area3D(), move.b->Area3D());

            bool regionReached = g->Count() <= targetAtlasSize;
            bool sizeThresholdReached = minNextArea > minChartArea;

            if (regionReached && sizeThresholdReached) {
                break;
            } else {
                moveCount++;
                bool accepted = ExecuteInitializedMove();
                if (accepted)
                    mergeCount++;
                if (mergeCount % 50 == 0)
                    LOG_VERBOSE << "Merged " << mergeCount << " regions...";
            }
        }

        ResetMoveState();

        LOG_INFO << "Clustering: Iteration "  << numIter << " took " << timer.TimeSinceLastCheck()
                 << " seconds. Attempts: " << moveCount << " Merges: " << mergeCount;

        numIter++;

    } while (mergeCount > 0);

    LOG_INFO << "Stopping after " << numIter << " passes and " << timer.TimeElapsed() << " seconds";

    LOG_DEBUG << "State dump";
    LOG_DEBUG << "Number of nodes: " << g->charts.size() << " (node list follows)";
    for (auto& entry : g->charts) {
        LOG_DEBUG << "Chart " << entry.first << " (FN = " << entry.second->FN() << ")";
    }
    LOG_DEBUG << "Number of moves: " << moves.size() << " (edge list follows)";
    for (const auto& entry : moves) {
        LOG_DEBUG << entry.first.a->id << " " << entry.first.b->id << " | w = " << entry.second;
    }

    LOG_DEBUG << " --- Execution stats ---";
    LOG_DEBUG << "\tAttempted moves: " << total_accepted + total_rejected;
    LOG_DEBUG << "\tAccepted moves: " << total_accepted;
    LOG_DEBUG << "\tRejected moves: " << total_rejected;
    LOG_DEBUG << "\tFailed feasibility checks: " << feasibility_failed;
    LOG_DEBUG << "\tFailed post checks: " << post_failed;
    LOG_DEBUG << "\tINIT time: " << init_t;
    LOG_DEBUG << "\tCHECK time: " << check_t;
    LOG_DEBUG << "\tMERGE time: " << merge_t;
    LOG_DEBUG << "\tOPT time: " << opt_t;

    LOG_DEBUG << "\t\tSHELL0 time: " << shell_0_t;
    LOG_DEBUG << "\t\tSHELL1 time: " << shell_1_t;
    LOG_DEBUG << "\t\tSHELL2 time: " << shell_2_t;
    LOG_DEBUG << "\t\tSHELL3 time: " << shell_3_t;
    LOG_DEBUG << "\t\tARAP time (total): " << total_arap_t;
    LOG_DEBUG << "\t\tARAP time (active): " << active_arap_t;

    LOG_DEBUG << "\tPOST time: " << post_t;

}

bool AtlasClustering::Valid(const WeightedClusteringMove& wm) const
{
    auto wmref = moves.find(wm.first);
    return wmref != moves.end() && wmref->second == wm.second;
}

bool AtlasClustering::InitializeNextMove()
{
    while (queue.size() > 0) {
        ensure_condition(state == MoveState::Uninitialized);
        WeightedClusteringMove wm = queue.top();
        if (Valid(wm)) {
            if (wm.second == Infinity()) {
                return false;
            } else {
                InitMove(wm.first);
                if (CheckFeasibility() == true) {
                    return true;
                } else {
                    RejectMove();
                    queue.pop();
                    InvalidateMove(wm.first);
               }
            }
        } else {
            queue.pop();
        }
    }
    return false;
}

bool AtlasClustering::ExecuteInitializedMove()
{
    ensure_condition(state == MoveState::Feasible);
    Merge();
    if (PostCheck()) {
        AcceptMove();
        return true;
    } else {
        RejectMove();
        return false;
    }
}
