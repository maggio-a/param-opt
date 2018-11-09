#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <cmath>
#include <Eigen/Core>

#ifndef M_PI
  #define M_PI 3.1415926535897932384626433832795
#endif
#ifndef M_PI_2
  #define M_PI_2 1.5707963267948966192313216916398
#endif

constexpr double Infinity()
{
    return std::numeric_limits<double>::infinity();
}

template <typename FaceType>
inline double EdgeLength(const FaceType& f, int i)
{
    return (f.cV0(i)->P() - f.cV1(i)->P()).Norm();
}

template <typename FaceType>
inline double EdgeLengthUV(const FaceType& f, int i)
{
    return (f.cWT(i).P() - f.cWT((i+1)%3).P()).Norm();
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
    //ensure_condition(v1.Norm() > 0 && v2.Norm() > 0);
    double v1n = v1.Norm();
    double v2n = v2.Norm();
    if (v1n == 0 || v2n == 0) {
        if (v1n == 0) v1n = 1e-6;
        if (v2n == 0) v2n = 1e-6;
    }
    double theta = VecAngle(v1, v2);
    if (!(theta > 0 && theta < M_PI)) {
        if (theta == 0) theta = 1e-3; // push theta to be very small
        else if (theta == M_PI) theta = M_PI - 1e-3; // almost flat
        else theta = M_PI / 2.0; // completely arbitrary value
    }
    w1[0] = v1n;
    w1[1] = 0;
    w2[0] = v2n * std::cos(theta);
    w2[1] = v2n * std::sin(theta);
}

template <typename PointType>
Eigen::Matrix2d ComputeTransformationMatrix(const PointType& x10, const PointType& x20, const PointType& u10, const PointType& u20)
{
    Eigen::Matrix2d f;
    Eigen::Matrix2d g;

    f(0, 0) = x10[0];
    f(1, 0) = x10[1];
    f(0, 1) = x20[0];
    f(1, 1) = x20[1];
    g(0, 0) = u10[0];
    g(1, 0) = u10[1];
    g(0, 1) = u20[0];
    g(1, 1) = u20[1];

    return g * f.inverse();
}

#endif // MATH_UTILS_H

