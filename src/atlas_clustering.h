#ifndef ATLAS_CLUSTERING_H
#define ATLAS_CLUSTERING_H

#include "mesh_graph.h"


struct ClusteringMove;


using WeightedClusteringMove = std::pair<ClusteringMove, double>;

using StitchOffsetVec = std::vector<std::pair<vcg::Point2d, double>>;

/* A move is just a pair of chart references. */
struct ClusteringMove {
    ChartHandle a;
    ChartHandle b;

    ClusteringMove(std::shared_ptr<FaceGroup> c1, std::shared_ptr<FaceGroup> c2) : a{c1}, b{c2}
    {
        ensure_condition(a != b);
        if (b->id > a->id)
            std::swap(a, b);
    }

    ClusteringMove() : a{nullptr}, b{nullptr}
    {
    }

    bool operator==(const ClusteringMove& other) const noexcept
    {
        return a == other.a && b == other.b;
    }

    bool IsNull() const
    {
        return (a == nullptr);
    }

    static ClusteringMove NullMove()
    {
        return ClusteringMove();
    }
};

struct ClusteringMoveHasher {
    std::size_t operator()(const ClusteringMove& m) const noexcept
    {
        std::size_t seed = 0;
        seed ^= std::hash<std::shared_ptr<FaceGroup>>()(m.a) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<std::shared_ptr<FaceGroup>>()(m.b) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

struct WeightedClusteringMoveComparator {
    int operator()(const WeightedClusteringMove& m1, const WeightedClusteringMove& m2) const
    {
        return m1.second > m2.second;
    }
};

/* A boundary chain holds a pair of lists of PosF objects These lists encode the
 * boundary shared between the two charts. The idea is that those sequences are
 * coherent wrt the topology of the boundary. Note that a chain may have 'jumps'
 * if the boundary is not continuous, and may also encode only a portion of the
 * whole boundary. */
struct BoundaryChain {
    ChartHandle a;
    ChartHandle b;
    std::vector<PosF> pos_a;
    std::vector<PosF> pos_b;
};

class AtlasClustering {

    /* A reference to the mesh graph */
    std::shared_ptr<MeshGraph> g;

    /* Priority queue for clustering moves. Elements in this queue may be no
     * longer valid due to merge operations (their weight might have changed, or
     * the graph edge itself may no longer exist if it referenced a region that
     * was absorbed by another), so a move in this queue is valid only if it
     * exists in the 'moves' map and the wieights match */
    std::priority_queue<WeightedClusteringMove, std::vector<WeightedClusteringMove>, WeightedClusteringMoveComparator> queue;

    /* This map stores for each move in the priority queue its 'real' weight.
     * Since obsolete entries are not removed from the priority queue, this map
     * is used to discard them when extracting the next best move. */
    std::unordered_map<ClusteringMove, double, ClusteringMoveHasher> moves;

    /* The boundary chains linked to the various moves */
    std::unordered_map<ClusteringMove, BoundaryChain, ClusteringMoveHasher> chains;

    /* The cumulative time spent solving ARAP systems */
    double total_arap_t;

    /* The cumulative time spent solving ARAP systems for ACCEPTED moves */
    double active_arap_t;


    //bool AddMove(ChartHandle c1, ChartHandle c2, bool replace);
    //bool AddMove(const WeightedClusteringMove& move, bool replace);

    /* Computes the full boundary chain and the cost of the move mv */
    bool ComputeWeightedMove(const ClusteringMove& mv, WeightedClusteringMove& wm, BoundaryChain& chain);

    /* Computes the cost of the move mv according to the given boundary chain */
    double ComputeCost(const ClusteringMove& mv, const BoundaryChain& chain);

    /* Adds a move to the collection of valid moves */
    bool AddMove(const WeightedClusteringMove& wm);

    /* Deletes a move from the collection of valid moves */
    void DeleteMove(const ClusteringMove& move);

    /* Sets the cost of m to Infinity() */
    void InvalidateMove(const ClusteringMove& m);



    /* The following functions control the merging process. The idea is that a
     * move is first set as the 'current move', and then the steps that from
     * it are performed in sequence. Logically a move is atomic, and its step
     * should be executed in isolation up to the point where the move is either
     * rejected or accepted. However, the single steps are separated so that a
     * GUI client can step into each section of the move.
     * I don't like this design... */

    enum MoveState {
        Uninitialized,
        Initialized,
        Feasible,
        Unfeasible,
        Merged,
        Acceptable,
        Unacceptable
    };

    /* The state of the current move */
    MoveState state;

    /* The current move. The invariant is that currentMove.IsNull() == false if
     * and only if since the last call to InitMove() neither AcceptMove() or
     * RejectMove() have been called. */
    ClusteringMove currentMove;

    /* The chart that remains after the merge */
    ChartHandle merge_a;

    /* The chart that is merged into target */
    ChartHandle merge_b;

    /* The texture id of merge_a */
    int texIdA;

    /* The texture id of merge_b */
    int texIdB;

    /* The list of stitched uv coordinates, and their tolerance offset */
    StitchOffsetVec offsetVec;

    /* Vectors of backup uv coordinates */
    std::vector<double> savedWtA;
    std::vector<double> savedWtB;

    double move_arap_t;
    double init_t;
    double check_t;
    double merge_t;
    double opt_t;
    double post_t;
    double shell_0_t;
    double shell_1_t;
    double shell_2_t;
    double shell_3_t;

    void ResetMoveState();

    /* Initializes the current move that the clustering algorithm is performing.
     * State change: Unitinialized -> Initialized */
    void InitMove(ClusteringMove move);

    /* Checks if the move is topologically feasible (whether or not the merge of
     * the two charts results in a topological disk
     * StateChange: Initialized -> Feasible | Unfeasible */
    bool CheckFeasibility();

    /* Merges the two charts into one.
     * State change: Feasible -> Merged */
    void Merge();

    /* Internal step, runs ARAP after merging and stitching */
    void OptimizeAfterMerge();

    /* Checks if the result of the merge is acceptable.
     * State change: Merged -> Acceptable | Unacceptable */
    bool PostCheck();

    /* Accepts the move. Note that a move can be forcefully accepted even if it
     * is unacceptable. This never happens in the 'headless' clustering
     * algorithm implementation, but can be done from a friend class such as the
     * GUI client.
     * State change: Acceptable | Unacceptable -> Uninitialized */
    void AcceptMove();

    /* Rejects the move, reverting any side effect that may have occurred while
     * performing the previous steps.
     * State change: Acceptable | Unacceptable | Unfeasible -> Uninitialized */
    void RejectMove();

    /* Returns the current move */
    const ClusteringMove& CurrentMove()
    {
        return currentMove;
    }

    /* Returns the current state */
    MoveState CurrentState()
    {
        return state;
    }

public:

    enum MergeStatus {
        Allowed = 0,
        ErrorDisconnected,
        ErrorUnfeasible
    };

    /* Tests if the merging together a collection of charts is feasible in the sense
     * that it yields a single connected component that does not require cutting to
     * be converted into a topological disk.
     * If the charts can be merged together, it returns a queue of charts that
     * defines a possible order of the moves that must be performed to merge the
     * whole collection. The moves can be performed by extracting the first chart,
     * and iteratively merging it with the following chart in the queue.
     * If the charts cannot be merged, it returns a proper error code and an empty
     * queue. */
    std::pair<MergeStatus,std::queue<ChartHandle>> MergeCheck(const std::vector<ChartHandle> charts);

    AtlasClustering(std::shared_ptr<MeshGraph> gptr);

    void Run(std::size_t targetAtlasSize, double smallRegionAreaThreshold);

    bool Valid(const WeightedClusteringMove& wm) const;

    /* Initializes the next topologically feasible move */
    bool InitializeNextMove();

    /* Executes the move that has been initialized and passed the feasibility
     * check. Returns a boolean indicating whether the move has been accepted */
    bool ExecuteInitializedMove();

};

#endif // ATLAS_CLUSTERING_H
