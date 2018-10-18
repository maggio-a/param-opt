#ifndef MESH_UTILS_H
#define MESH_UTILS_H

#include "mesh.h"

#include <vector>

/* A set of functions to introduce cuts in shell meshes. Edges that are texture
 * seams in the input mesh are marked as 'FAUX' edges, and the shell meshes can
 * be cut along such edges to reduce distortion in the new parameterization. */

using PosF = face::Pos<MeshFace>;

struct PosNode {
    PosF pos;
    double distance;

    PosNode() = default;
    PosNode(PosF p, double d) : pos{p}, distance{d} {}

    bool operator<(const PosNode& other) const { return distance < other.distance; }
};

/* Closes the holes of a mesh using the vcg MinimumWeightEar strategy.
 * Hole-filling faces are marked by setting the holeFilling field to true */
void CloseMeshHoles(Mesh& mesh);

/* Marks initial texture seam edges as faux */
void MarkInitialSeamsAsFaux(Mesh& shell, Mesh& baseMesh);

/* Computes a fan of faux pos objects. This function is used when walking along
 * a faux edge path to retrieve the allowed next steps. */
std::vector<PosF> GetFauxPosFan(PosF& startPos);

/* Computes vertices distances to the boundary (restricted only to vertices that
 * lie on a seam) as vertex quality. The distance of vertices that are not
 * incident to seam edges is INFINITY. Assumes the seams are marked as faux.
 * NOTE (FIXME?) as of now distances are computed according to the target shape
 * attribute (does it make sense?) */
void ComputeDistanceFromBorderOnSeams(Mesh& m);

/* Selects the shortest path from a starting (seam) pos to the boundary, along
 * faux edges. This function allows to select the path for the subsequent cut.
 * Assumes the seams are marked as faux
 * Returns a Pos pointing to the reached boundary vertex */
PosF SelectShortestSeamPathToBoundary(Mesh& m, const PosF& pos);

/* Selects the shortest path from a starting (seam) pos to a peak, along faux
 * (seam) edges. A peak vertex is a vertex that with maximal finite distance
 * from the boundary */
void SelectShortestSeamPathToPeak(Mesh &m, const PosF& pos);

/* If the path touches any hole-filling area, it leaves selected only the path
 * that joins the boundary with the nearest such area.
 * Returns whether or not the path reached a fill area */
bool RectifyCut(Mesh& shell, PosF boundaryPos);

/* To be called after a shell has been cut, removes any obsolete faux flag and
 * deletes any hole-filling region that is adjacent to the boundary after a cut
 * has been applied
 * NOTE WARNING FIXME it relies on the edge selection flags, so they must NOT be
 * cleared before invoking this function */
void CleanupShell(Mesh& shell);

/* Copies a shell object with all its attributes */
void CopyShell(Mesh& shell, Mesh& out);

#endif // MESH_UTILS_H

