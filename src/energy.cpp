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
    : m{mesh}
{
    ensure_condition(HasTargetShapeAttribute(mesh));
    targetShape = GetTargetShapeAttribute(mesh);
}

Energy::~Energy()
{
    // empty
}

void Energy::UpdateCache()
{
    numMeshFaces = 0;
    numHoleFillingFaces = 0;
    numScaffoldFaces = 0;
    for (auto& f : m.face) {
        if (f.IsMesh())
            numMeshFaces++;
        else if (f.IsHoleFilling())
            numHoleFillingFaces++;
        else if (f.IsScaffold())
            numScaffoldFaces++;
        else
            ensure_condition(0);
    }
    ensure_condition(numMeshFaces + numHoleFillingFaces + numScaffoldFaces == (int) m.face.size());

    surfaceArea = 0;
    holeFillingArea = 0;
    for (auto& f : m.face) {
        double area = FaceArea(&f);
        surfaceArea += area;
        if (!f.IsMesh())
            holeFillingArea += area;
    }

    // This must be done after everything else has been updated, so it is 'safe'
    // to call E_IgnoreMarkedFaces() from the subclasses
    if (numScaffoldFaces > 0)
        scaffoldRegularizationTerm = (1.0 / 100.0) * (E_IgnoreMarkedFaces() / numScaffoldFaces);
    else
        scaffoldRegularizationTerm = 0;
}


double Energy::E()
{
    double totalEnergy = 0;
    for (auto& f : m.face) {
        double e = E(f);
        if (f.IsScaffold())
            totalEnergy += scaffoldRegularizationTerm * ((e / FaceArea(&f)) - NormalizedMinValue());
        else
            totalEnergy += e;
    }
    return totalEnergy;
}

double Energy::E_IgnoreMarkedFaces()
{
    double totalEnergy = 0;
    for (auto& f : m.face)
        if (f.IsMesh())
            totalEnergy += E(f);
    return totalEnergy;
}

MatrixXd Energy::Grad()
{
    MatrixXd g = MatrixXd::Zero(m.VN(), 2);
    for (auto& f : m.face) {
        Eigen::Vector2d gf[3];
        Grad(tri::Index(m, f), gf[0], gf[1], gf[2]);
        if (f.IsScaffold()) {
            double fa = FaceArea(&f);
            gf[0] *= (scaffoldRegularizationTerm / fa);
            gf[1] *= (scaffoldRegularizationTerm / fa);
            gf[2] *= (scaffoldRegularizationTerm / fa);
        }
        g.row(Index(m, f.V(0))) += gf[0];
        g.row(Index(m, f.V(1))) += gf[1];
        g.row(Index(m, f.V(2))) += gf[2];
    }
    return g;
}

void Energy::MapToFaceQuality(bool normalized)
{
    for (auto& f : m.face) {
        f.Q() = E(f);
        if (f.IsScaffold())
            f.Q() /= (scaffoldRegularizationTerm / FaceArea(&f));
        else if (normalized)
            f.Q() /= FaceArea(&f);
    }
}

double Energy::FaceArea(Mesh::ConstFacePointer fp)
{
    return (((P(fp, 1) - P(fp, 0)) ^ (P(fp, 2) - P(fp, 0))).Norm() / 2.0);
}

double Energy::SurfaceArea()
{
    return surfaceArea;
}

double Energy::SurfaceAreaNotHoleFilling()
{
    return surfaceArea - holeFillingArea;
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

SymmetricDirichletEnergy::SymmetricDirichletEnergy(Mesh& mesh)
    : Energy{mesh}, data{m.face}
{
    UpdateCache();
}

double SymmetricDirichletEnergy::NormalizedMinValue()
{
    return 4.0;
}

double SymmetricDirichletEnergy::E(const Mesh::FaceType& f)
{
    double areaUV = 0.5 * ((u1 - u0) ^ (u2 - u0));

    if (areaUV <= 0) {
        return Infinity();
    } else {
        double area3D = data[f][3];
        double o[3] = { // (opposite edge)^2
            (   u1-u2).SquaredNorm(),
            (u0   -u2).SquaredNorm(),
            (u0-u1   ).SquaredNorm(),
        };
        double e_d = 0.5 * (data[f][0] * o[0] + data[f][1] * o[1] + data[f][2] * o[2]);
        double energy = (1 + (area3D*area3D)/(areaUV*areaUV)) * (e_d);
        return energy;
    }
}

void SymmetricDirichletEnergy::Grad(int faceIndex, Eigen::Vector2d& g0, Eigen::Vector2d& g1, Eigen::Vector2d& g2)
{
        const MeshFace& f = m.face[faceIndex];

        double area3D = data[f][3];
        double areaUV = 0.5 * ((u1 - u0) ^ (u2 - u0));

        double angleTermLeft = (u2-u0).SquaredNorm() * (P(&f, 1) - P(&f, 0)).SquaredNorm() + (u1-u0).SquaredNorm() * (P(&f, 2) - P(&f, 0)).SquaredNorm();
        angleTermLeft /= 4.0 * area3D;

        double angleTermRight = ((u2-u0) * (u1-u0)) * ((P(&f, 2) - P(&f, 0)) * (P(&f, 1) - P(&f, 0)));
        angleTermRight /= 2.0 * area3D;

        double angleTerm = angleTermLeft - angleTermRight;

        double areaTermRatio = - (area3D*area3D) / std::pow(areaUV, 3);
        double e_a = (1 + (area3D*area3D) / (areaUV*areaUV));

        Eigen::Vector2d g[3];

        for (int i = 0; i < 3; ++i) {
            int j = (i+1)%3;
            int k = (i+2)%3;

            double gu_area = areaTermRatio * ( uj.Y() - uk.Y());
            double gv_area = areaTermRatio * (-uj.X() + uk.X());

            double gu_angle = (data[f][k]*(ui.X() - uj.X()) + data[f][j]*(ui.X() - uk.X()));
            double gv_angle = (data[f][k]*(ui.Y() - uj.Y()) + data[f][j]*(ui.Y() - uk.Y()));

            double gu = gu_area * angleTerm + e_a * gu_angle;
            double gv = gv_area * angleTerm + e_a * gv_angle;

            ensure_condition(std::isnan(gu) == false && std::isnan(gv) == false);

            g[i] = Eigen::Vector2d{gu, gv};
        }

        g0 = g[0]; g1 = g[1]; g2 = g[2];
}

void SymmetricDirichletEnergy::UpdateCache()
{
    data.UpdateSize();
    for (auto&f : m.face) {
        data[f][3] = (((P(&f, 1) - P(&f, 0)) ^ (P(&f, 2) - P(&f, 0))).Norm() / 2.0);
        ensure_condition(!std::isnan(data[f][3]));
        ensure_condition(std::isfinite(data[f][3]));
        ensure_condition(data[f][3] > 0);
    }
    for (auto& f : m.face) {
        for (int i=0; i<3; i++) {
            Point3d a = P1(&f, i) - P0(&f, i);
            Point3d c = P2(&f, i) - P0(&f, i);
            data[f][i] = VecCotg(a, c);
            ensure_condition(!std::isnan(data[f][i]));
            ensure_condition(std::isfinite(data[f][i]));
        }
    }
    Energy::UpdateCache();
}

#undef u0
#undef u1
#undef u2
#undef ui
#undef uj
#undef uk
