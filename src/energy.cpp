#include "energy.h"
#include "math_utils.h"
#include "metric.h"
#include "mesh_attribute.h"

#include <Eigen/Core>


#define u0 (f.cV(0)->T().P())
#define u1 (f.cV(1)->T().P())
#define u2 (f.cV(2)->T().P())
#define ui (f.cV(i)->T().P())
#define uj (f.cV(j)->T().P())
#define uk (f.cV(k)->T().P())


using Eigen::MatrixXd;

/* Energy interface implementation
 * =============================== */

Energy::Energy(Mesh& mesh)
    : m{mesh},
      surfaceArea{0.0},
      holeFillingArea{0.0}
{
    assert(HasTargetShapeAttribute(mesh));
    targetShape = GetTargetShapeAttribute(mesh);

    for (auto& f : m.face) {
        double area = (((P(&f, 1) - P(&f, 0)) ^ (P(&f, 2) - P(&f, 0))).Norm() / 2.0);
        surfaceArea += area;
        if (f.holeFilling) holeFillingArea += area;
    }
}

Energy::~Energy()
{
    // empty
}

double Energy::E()
{
    double e = 0;
    for (auto& f : m.face) e += E(f);
    return e;
}

double Energy::E_IgnoreMarkedFaces(bool normalized)
{
    double e = 0;
    for (auto& f : m.face)
        if (f.holeFilling == false)
            e += E(f);
    return e / ((normalized) ? (surfaceArea - holeFillingArea) : 1.0);
}

void Energy::UpdateCache()
{
    // empty
}

void Energy::MapToFaceQuality(bool normalized)
{
    for (auto& f : m.face) {
        f.Q() = E(f, normalized);
    }
}

Mesh::CoordType Energy::P(Mesh::ConstFacePointer fp, int i)
{
    assert(i >= 0 && i <= 2);
    return targetShape[fp].P[i];
}

double Energy::FaceArea(Mesh::ConstFacePointer fp)
{
    return (((P(fp, 1) - P(fp, 0)) ^ (P(fp, 2) - P(fp, 0))).Norm() / 2.0);
}

double Energy::SurfaceArea()
{
    return surfaceArea;
}

double Energy::ParameterArea()
{
    double parameterArea = 0.0;
    for (auto& f : m.face) {
        parameterArea += 0.5 * ((u1 - u0) ^ (u2 - u0));
    }
    return parameterArea;
}

void Energy::CorrectScale()
{
    double scale = std::sqrt(SurfaceArea() / ParameterArea());
    for (auto& v : m.vert) {
        v.T().P() = scale * v.T().P();
    }
}


/* Symmetric Dirichlet energy implementation
 * ========================================= */

SymmetricDirichlet::SymmetricDirichlet(Mesh& mesh)
    : Energy{mesh}, data{m.face}
{
    UpdateCache();
}

double SymmetricDirichlet::E(const Mesh::FaceType& f, bool normalized)
{
    double o[3] = { // (opposite edge)^2
        (   u1-u2).SquaredNorm(),
        (u0   -u2).SquaredNorm(),
        (u0-u1   ).SquaredNorm(),
    };

    double area3D = data[f][3];
    double areaUV = 0.5 * ((u1 - u0) ^ (u2 - u0));
    double e_d = 0.5 * (data[f][0] * o[0] + data[f][1] * o[1] + data[f][2] * o[2]);

    double energy = (1 + (area3D*area3D)/(areaUV*areaUV)) * (e_d);

    if (normalized) energy /= area3D;

    double attenuation = 1.0;
    if (f.holeFilling) attenuation = 0.0001;
    return attenuation * energy;
}

MatrixXd SymmetricDirichlet::Grad()
{
    MatrixXd g = MatrixXd::Zero(m.VN(), 2);

    for(auto& f : m.face) {
        double area3D = data[f][3];
        double areaUV = 0.5 * ((u1 - u0) ^ (u2 - u0));

        double angleTermLeft = (u2-u0).SquaredNorm() * (P(&f, 1) - P(&f, 0)).SquaredNorm() + (u1-u0).SquaredNorm() * (P(&f, 2) - P(&f, 0)).SquaredNorm();
        angleTermLeft /= 4.0 * area3D;

        double angleTermRight = ((u2-u0) * (u1-u0)) * ((P(&f, 2) - P(&f, 0)) * (P(&f, 1) - P(&f, 0)));
        angleTermRight /= 2.0 * area3D;

        double angleTerm = angleTermLeft - angleTermRight;

        double areaTermRatio = - (area3D*area3D) / std::pow(areaUV, 3);
        double e_a = (1 + (area3D*area3D) / (areaUV*areaUV));

        double attenuation = 1.0;
        if (f.holeFilling) attenuation = 0.0001;

        for (int i = 0; i < 3; ++i) {
            int j = (i+1)%3;
            int k = (i+2)%3;

            double gu_area = areaTermRatio * ( uj.Y() - uk.Y());
            double gv_area = areaTermRatio * (-uj.X() + uk.X());

            double gu_angle = (data[f][k]*(ui.X() - uj.X()) + data[f][j]*(ui.X() - uk.X()));
            double gv_angle = (data[f][k]*(ui.Y() - uj.Y()) + data[f][j]*(ui.Y() - uk.Y()));

            double gu = gu_area * angleTerm + e_a * gu_angle;
            double gv = gv_area * angleTerm + e_a * gv_angle;

            assert(std::isnan(gu) == false && std::isnan(gv) == false);

            g.row(Index(m, f.V(i))) += attenuation * Eigen::Vector2d{gu, gv};
        }
    }

    return g;
}

void SymmetricDirichlet::UpdateCache()
{
    data.UpdateSize();
    for (auto&f : m.face) {
        data[f][3] = (((P(&f, 1) - P(&f, 0)) ^ (P(&f, 2) - P(&f, 0))).Norm() / 2.0);
        assert(data[f][3] > 0);
    }
    for (auto& f : m.face) {
        for (int i=0; i<3; i++) {
            Point3d a = P1(&f, i) - P0(&f, i);
            Point3d c = P2(&f, i) - P0(&f, i);
            data[f][i] = VecCotg(a, c);
        }
    }
}

#undef u0
#undef u1
#undef u2
#undef ui
#undef uj
#undef uk
