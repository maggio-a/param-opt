#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <cmath>

template <typename FaceType>
inline double EdgeLength(const FaceType& f, int i)
{
    return (f.cV0(i)->P() - f.cV1(i)->P()).Norm();
}


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

/* Given two vectors, it transforms them to the local 2d-frame of the plane they span */
template <typename PointType, typename PointTypeOut>
void LocalIsometry(const PointType& v1, const PointType& v2, PointTypeOut& w1, PointTypeOut& w2)
{
    double theta = VecAngle(v1, v2);
    double v2Norm = v2.Norm();
    w1[0] = v1.Norm();
    w1[1] = 0;
    w2[0] = v2Norm * std::cos(theta);
    w2[1] = v2Norm * std::sin(theta);
}

#endif // MATH_UTILS_H

