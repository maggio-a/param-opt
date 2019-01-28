#ifndef CONES_H
#define CONES_H

class Mesh;

/* returns the index of the vertex that is subject to the largest area stretch
 * under a conformal transformation */
int FindApproximateCone(Mesh& m, bool selectedFlag);

#endif // CONES_H
