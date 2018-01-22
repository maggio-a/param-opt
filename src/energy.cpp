#include "energy.h"
#include "math_utils.h"

#include <Eigen/Core>


#define p0 (f.cV(0)->P())
#define p1 (f.cV(1)->P())
#define p2 (f.cV(2)->P())
#define u0 (f.cV(0)->T().P())
#define u1 (f.cV(1)->T().P())
#define u2 (f.cV(2)->T().P())
#define ui (f.cV(i)->T().P())
#define uj (f.cV(j)->T().P())
#define uk (f.cV(k)->T().P())


using Eigen::MatrixXd;


Energy::Energy(Mesh& mesh, Geometry geometryMode) : m{mesh}, mode{geometryMode}, surfaceArea{0.0}
{
    tcsattr = tri::Allocator<Mesh>::FindPerFaceAttribute<TexCoordStorage>(m, "WedgeTexCoordStorage");
    assert(tri::Allocator<Mesh>::IsValidHandle<TexCoordStorage>(m, tcsattr));
    for (auto& f : m.face) {
        surfaceArea += (((P(&f, 1) - P(&f, 0)) ^ (P(&f, 2) - P(&f, 0))).Norm() / 2.0);
    }
}

Mesh::CoordType Energy::P(Mesh::ConstFacePointer fp, int i) {
    assert(i>=0 && i<3);
    if (mode == Geometry::Model) return fp->cP(i);
    else if (mode == Geometry::Texture) return Mesh::CoordType(tcsattr[fp].tc[i].U(), tcsattr[fp].tc[i].V(), 0);
    else {
        assert(0 && "Energy::P()");
        return {0, 0, 0};
    }
}

Mesh::CoordType Energy::P0(Mesh::ConstFacePointer fp, int i) { return P(fp, i); }
Mesh::CoordType Energy::P1(Mesh::ConstFacePointer fp, int i) { return P(fp, (i+1)%3); }
Mesh::CoordType Energy::P2(Mesh::ConstFacePointer fp, int i) { return P(fp, (i+2)%3); }

double Energy::SurfaceArea() {
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


// Symmetric Dirichlet energy implementation
// =========================================

SymmetricDirichlet::SymmetricDirichlet(Mesh& mesh, Geometry geometryMode = Geometry::Model) : Energy{mesh, geometryMode}, data{m.face}
{
    for (auto&f : m.face) {
        data[f][3] = (((P(&f, 1) - P(&f, 0)) ^ (P(&f, 2) - P(&f, 0))).Norm() / 2.0f);
    }

    for (auto& f : m.face) {
        for (int i=0; i<3; i++) {
            // Numerically stable (?) cotangents
            Point3d a = P1(&f, i) - P0(&f, i);
            //Point3d b = P2(&f, i) - P1(&f, i);
            Point3d c = P2(&f, i) - P0(&f, i);
            //double cotg = (a.SquaredNorm() + c.SquaredNorm() - b.SquaredNorm())/(2.0*data[f][3])/2.0;
            //data[f][i] = cotg;
            data[f][i] = VecCotg(a, c);
        }
    }
}

double SymmetricDirichlet::E()
{
    double e = 0;
    for (auto& f : m.face) e += E(f);
    return e;
}

double SymmetricDirichlet::E(const Mesh::FaceType& f)
{
    double o[3] = { // (opposite edge)^2
        (   u1-u2).SquaredNorm(),
        (u0   -u2).SquaredNorm(),
        (u0-u1   ).SquaredNorm(),
    };

    double area3D = data[f][3];
    double areaUV = 0.5 * ((u1 - u0) ^ (u2 - u0));

    double e_d = 0.5 * (data[f][0] * o[0] + data[f][1] * o[1] + data[f][2] * o[2]);
    //double e_d = (data[f][0] * o[0] + data[f][1] * o[1] + data[f][2] * o[2]) / (2.0 * area3D);
    //double a = std::atan(1.0 / data[f][0]);
    //double b = std::atan(1.0 / data[f][1]);
    //double c = std::atan(1.0 / data[f][2]);

    /*double angle[3];
    Point2d coord[3];
    for (int i=0; i<3; i++) {
        coord[i] = Point2d(P(&f, i)[0], P(&f, i)[1]);
        angle[i] = std::max(vcg::Angle(P1(&f, i) - P0(&f, i), P2(&f, i) - P0(&f, i)), float(1e-8));
    }*/

    //double test = (coord[1] - coord[0]) ^ (coord[2] - coord[0]);

    //Point3f pa = P0(&f);
    //Point3f pb = P1(&f);
    //Point3f pc = P2(&f);


    //double energy = area3D * (1 + (area3D*area3D)/(areaUV*areaUV)) * (e_d);
    //double eee = (1 + (area3D*area3D)/(areaUV*areaUV));
    double energy = (1 + (area3D*area3D)/(areaUV*areaUV)) * (e_d);
    //assert (e_d >= 0);
    return energy;
}

MatrixXd SymmetricDirichlet::Grad()
{
    MatrixXd g = MatrixXd::Zero(m.VN(), 2);

    for(auto& f : m.face) {
        double o[3] = { // (opposite edge)^2
            (   u1-u2).SquaredNorm(),
            (u0   -u2).SquaredNorm(),
            (u0-u1   ).SquaredNorm(),
        };

        double area3D = data[f][3];
        double areaUV = 0.5 * ((u1 - u0) ^ (u2 - u0));

        //double e_d = 0.5 * (data[f][0] * o[0] + data[f][1] * o[1] + data[f][2] * o[2]);
        double angleTermLeft = (u2-u0).SquaredNorm() * (P(&f, 1) - P(&f, 0)).SquaredNorm() + (u1-u0).SquaredNorm() * (P(&f, 2) - P(&f, 0)).SquaredNorm();
        angleTermLeft /= 4.0 * area3D;

        double angleTermRight = ((u2-u0) * (u1-u0)) * ((P(&f, 2) - P(&f, 0)) * (P(&f, 1) - P(&f, 0)));
        angleTermRight /= 2.0 * area3D;

        double angleTerm = angleTermLeft - angleTermRight;

        double areaTermRatio = - (area3D*area3D) / std::pow(areaUV, 3);
        double e_a = (1 + (area3D*area3D) / (areaUV*areaUV));

        //assert(areaUV > 0);

        for (int i = 0; i < 3; ++i) {
            int j = (i+1)%3;
            int k = (i+2)%3;

            double gu_area = areaTermRatio * ( uj.Y() - uk.Y());
            double gv_area = areaTermRatio * (-uj.X() + uk.X());

            double gu_angle = (data[f][k]*(ui.X() - uj.X()) + data[f][j]*(ui.X() - uk.X()));
            double gv_angle = (data[f][k]*(ui.Y() - uj.Y()) + data[f][j]*(ui.Y() - uk.Y()));

            //double gu = gu_area * e_d + e_a * gu_angle;
            double gu = gu_area * angleTerm + e_a * gu_angle;
            //double gv = gv_area * e_d + e_a * gv_angle;
            double gv = gv_area * angleTerm + e_a * gv_angle;

            g.row(Index(m, f.V(i))) += Eigen::Vector2d{gu, gv};
        }
    }

    return g;
}

#undef p0
#undef p1
#undef p2
#undef u0
#undef u1
#undef u2
#undef ui
#undef uj
#undef uk
