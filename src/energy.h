#ifndef ENERGY_H
#define ENERGY_H

#include "mesh.h"
#include "uv.h"

#include <Eigen/Core>


class Energy {

    friend class DescentMethod;
    friend class SLIM;

    Mesh::PerFaceAttributeHandle<TexCoordStorage> tcsattr;

public:

    enum Geometry { Model, Texture };

protected:

    Mesh& m;
    Geometry mode;
    double surfaceArea;

public:

    Energy(Mesh& mesh, Geometry geometryMode);
    virtual ~Energy() {}

    virtual double E() = 0;
    virtual double E(const Mesh::FaceType& f) = 0;
    virtual Eigen::MatrixXd Grad() = 0;

    double SurfaceArea();
    double ParameterArea();
    void CorrectScale();

    Mesh::CoordType P(Mesh::ConstFacePointer fp, int i);
    Mesh::CoordType P0(Mesh::ConstFacePointer fp, int i = 0);
    Mesh::CoordType P1(Mesh::ConstFacePointer fp, int i = 0);
    Mesh::CoordType P2(Mesh::ConstFacePointer fp, int i = 0);

};


class SymmetricDirichlet : public Energy {

    SimpleTempData<Mesh::FaceContainer, Point4d> data; // cotangents and surface area

public:

    SymmetricDirichlet(Mesh& mesh, Geometry geometryMode);

    double E();
    double E(const Mesh::FaceType& f);

    Eigen::MatrixXd Grad();
};


#endif // ENERGY_H

