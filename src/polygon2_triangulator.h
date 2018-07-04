#ifndef POLYGON2_TRIANGULATOR_H
#define POLYGON2_TRIANGULATOR_H

//#include <mapbox/earcut.hpp>

#include <array>
#include <vector>
#include <cassert>

#define ANSI_DECLARATORS
#define REAL double
#define VOID void
extern "C" {
    #include <triangle.h>
}
#undef ANSI_DECLARATORS
#undef REAL
#undef VOID

using Poly2Point = std::array<double,2>;
using Polyline2 = std::vector<std::array<double,2>>;
using Poly2 = std::vector<Polyline2>;

inline std::array<double,2> Poly2PointLerp(const Poly2Point& p1, const Poly2Point& p2, double t)
{
    assert(t >= 0 && t <= 1);
    return { ((1.0 - t) * p1[0] + t * p2[0]), ((1.0 - t) * p1[1] + t * p2[1]) };
}

/*
inline Polyline2 BuildPolyline2(const std::vector<std::size_t> &vfi, const std::vector<int> &vvi, const Mesh& shell2D)
{
    Polyline2 polyline;
    polyline.reserve(vfi.size());
    for (std::size_t i = 0; i < vfi.size(); ++i) {
        const vcg::Point3d& p = shell2D.face[vfi[i]].cP(vvi[i]);
        assert(p.Z() == 0.0);
        polyline.push_back({p.X(), p.Y()});
    }
    return polyline;
}
*/

/*
inline double LenPolyline2(const std::vector<std::size_t> &vfi, const std::vector<int> &vvi, const Mesh& shell2D)
{
    std::vector<vcg::Point3d> polyline;
    double len = 0;
    for (std::size_t i = 0; i < vfi.size(); ++i) {
        polyline.push_back(shell2D.face[vfi[i]].cP(vvi[i]));
    }
    for (std::size_t i = 1; i < polyline.size(); ++i) {
        len += (polyline[i] - polyline[i-1]).Norm();
    }
    len += (polyline.front() - polyline.back()).Norm();
    return len;
}
*/

inline void Triangulate(const Poly2& poly, const std::vector<double>& holes, Polyline2& pointsOut, std::vector<unsigned>& indices)
{
    assert(holes.size() % 2 == 0);
    pointsOut.clear();
    indices.clear();
    struct triangulateio in = {};  // zero
    struct triangulateio out = {}; // zero

    //memset(&in, 0, sizeof(struct triangulateio));
    //memset(&out, 0, sizeof(struct triangulateio));

    // Merge all the Polyline2 vertices in a global list
    std::vector<double> points;
    for (const auto& polyline : poly) {
        for (const auto& p : polyline) {
            points.push_back(p[0]);
            points.push_back(p[1]);
        }
    }
    // Treat each Polyline2 as a list of segments
    std::vector<int> segs;
    for (const auto& polyline : poly) {
        int segbase = segs.size() / 2; // index of the first vertex of the segment in the global list
        for (unsigned k = 0; k < polyline.size(); k++) {
            segs.push_back(segbase + k);
            segs.push_back(segbase + ((k+1) % polyline.size()));
        }
    }

    std::vector<double> holescopy = holes;

    // Fill the input data structure
    in.pointlist = &points[0];
    in.numberofpoints = points.size() / 2;
    in.segmentlist = &segs[0];
    in.numberofsegments = segs.size() / 2;
    in.holelist = &holescopy[0];
    in.numberofholes = holes.size() / 2;

    char args[] = "zpqYYQ";
    triangulate(args, &in, &out, NULL);

    pointsOut.reserve(out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; i++)
        pointsOut.push_back({*(out.pointlist + (2 * i)), *(out.pointlist + (2 * i) + 1)});

    indices.reserve(out.numberoftriangles);
    for (int i = 0; i < out.numberoftriangles; ++i) {
        indices.push_back(*(out.trianglelist + (3 * i)));
        indices.push_back(*(out.trianglelist + (3 * i) + 1));
        indices.push_back(*(out.trianglelist + (3 * i) + 2));
    }

    // cleanup
    free(out.pointlist);
    free(out.pointattributelist);
    free(out.pointmarkerlist);
    free(out.trianglelist);
    free(out.triangleattributelist);
    free(out.neighborlist);
    free(out.segmentlist);
    free(out.segmentmarkerlist);
    free(out.edgelist);
    free(out.edgemarkerlist);
}


#endif // POLYGON2_TRIANGULATOR_H

