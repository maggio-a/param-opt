#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <cmath>

/* Computes the angle between u and v */
template <typename PointType>
double VecAngle(const PointType& u, const PointType& v)
{
   typename PointType::ScalarType nu = u.Norm();
   typename PointType::ScalarType nv = v.Norm();

   double n = (u*nv - v*nu).Norm();
   double d = (u*nv + v*nu).Norm();

   return 2.0 * std::atan(n/d);
}

/* Computes the cotangent of the angle between u and v */
template <typename PointType>
double VecCotg(const PointType& u, const PointType& v)
{
    const PointType w = v - u;
    double dblArea = (u ^ v).Norm();
    return (u.SquaredNorm() + v.SquaredNorm() - w.SquaredNorm()) / (2.0 * dblArea);
}

#endif // MATH_UTILS_H

