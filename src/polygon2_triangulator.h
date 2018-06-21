#ifndef POLYGON2_TRIANGULATOR_H
#define POLYGON2_TRIANGULATOR_H

#include <mapbox/earcut.hpp>

using Poly2Point = std::array<double,2>;
using Polyline2 = std::vector<std::array<double,2>>;
using Poly2 = std::vector<Polyline2>;

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


#endif // POLYGON2_TRIANGULATOR_H

