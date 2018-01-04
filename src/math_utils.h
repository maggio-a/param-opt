#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <cmath>

template <typename PointType>
static double VecAngle(const PointType& a, const PointType& b)
{
   typename PointType::ScalarType na = a.Norm();
   typename PointType::ScalarType nb = b.Norm();

   double n = (a*nb - b*na).Norm();
   double d = (a*nb + b*na).Norm();

   return 2.0 * std::atan(n/d);
}

#endif // MATH_UTILS_H

