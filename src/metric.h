#ifndef METRIC_H
#define METRIC_H

#include "mesh.h"
#include "vertex_position.h"
#include "math_utils.h"

enum ParameterizationGeometry { Model, Texture };

class DistortionMetric {

public:

    enum Type { Area, Angle };

    static void ComputeAreaScale(Mesh& m, double& areaScale, ParameterizationGeometry geometry)
    {
        double sumArea3D = 0;
        double sumAreaUV = 0;
        for (const auto& f : m.face) {
            double a3D = Area3D(m, f, geometry);
            double aUV = AreaUV(f);
            if (a3D > 0 && aUV > 0) {
                sumArea3D += a3D;
                sumAreaUV += aUV;
            }
        }
        areaScale = sumArea3D / sumAreaUV;
    }

    static double AreaUV(const Mesh::FaceType& f)
    {
        Point2d u0 = f.cWT(0).P();
        Point2d u1 = f.cWT(1).P();
        Point2d u2 = f.cWT(2).P();
        return ((u1 - u0) ^ (u2 - u0)) / 2.0;
    }

    static double AngleUV(const Mesh::FaceType& f, int i)
    {
        Point2d u0 = f.cWT(i).P();
        Point2d u1 = f.cWT((i+1)%3).P();
        Point2d u2 = f.cWT((i+2)%3).P();
        return VecAngle(u1 - u0, u2 - u0);
    }

    static double Angle3D(Mesh& m, const Mesh::FaceType& f, int i, ParameterizationGeometry geometry)
    {
        int j = (i+1)%3;
        int k = (i+2)%3;
        if (geometry == Model) {
            return VecAngle(f.cP(j) - f.cP(i), f.cP(k) - f.cP(i));
        } else {
            WedgeTexCoordAttributePosition<Mesh> vpos{m, "WedgeTexCoordStorage"};
            return VecAngle<Point3d>(vpos(&f, j) - vpos(&f, i), vpos(&f, k) - vpos(&f, i));
        }
    }

    static double Area3D(const Mesh::FaceType& f)
    {
        return ((f.cP(1) - f.cP(0)) ^ (f.cP(2) - f.cP(0))).Norm() / 2.0;
    }

    static double Area3D(Mesh& m, const Mesh::FaceType& f, ParameterizationGeometry geometry)
    {
        if (geometry == Model) {
            return ((f.cP(1) - f.cP(0)) ^ (f.cP(2) - f.cP(0))).Norm() / 2.0;
        } else {
            WedgeTexCoordAttributePosition<Mesh> vpos{m, "WedgeTexCoordStorage"};
            return ((vpos(&f, 1) - vpos(&f, 0)) ^ (vpos(&f, 2) - vpos(&f, 0))).Norm() / 2.0;
        }
    }

    static double AreaDistortion(Mesh& m, const Mesh::FaceType& f, double areaScale, ParameterizationGeometry geometry)
    {
        double parameterArea = std::abs(AreaUV(f) * areaScale);
        double modelArea = std::abs(Area3D(m, f, geometry));
        return (parameterArea - modelArea) / modelArea;
    }

    static double AngleDistortion(Mesh& m, const Mesh::FaceType& f, ParameterizationGeometry geometry)
    {
        double d = 0;
        for (int i = 0; i < 3; ++i) {
            double parameterAngle = AngleUV(f, i);
            double faceAngle = Angle3D(m, f, i, geometry);
            d += std::abs(parameterAngle - faceAngle);
        }
        return d;
    }
};


#endif

