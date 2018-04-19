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

    enum Geometry { Model, Texture, Mixed };

protected:

    Mesh& m;
    const Geometry mode;
    double surfaceArea;
    double holeFillingArea;
    SimpleTempData<Mesh::FaceContainer, Point3d> P0_; // cotangents and surface area
    SimpleTempData<Mesh::FaceContainer, Point3d> P1_; // cotangents and surface area
    SimpleTempData<Mesh::FaceContainer, Point3d> P2_; // cotangents and surface area


public:

    Energy(Mesh& mesh, Geometry geometryMode);
    virtual ~Energy() {}

    double E();
    double E_IgnoreMarkedFaces(bool normalized = false);
    virtual double E(const Mesh::FaceType& f) = 0;
    virtual Eigen::MatrixXd Grad() = 0;

    double FaceArea(Mesh::ConstFacePointer fp);
    double SurfaceArea();
    double ParameterArea();

     /*
      * Hack to speed up the convergence of iterative methods
      *
      * Scales the parameterization so that its area matches the total area of the model (either the actual 3D area or
      * the area of the original parameterized faces, according to the geometry mode).
      * This can be used by iterative methods on energies that penalize area mismatches, in order to scale the initial solution
      * and reduce the scaling process that happens during the iterations (thereby lowering the number of iterations required
      * to converge)
      * */
    void CorrectScale();

    Mesh::CoordType P(Mesh::ConstFacePointer fp, int i);
    Mesh::CoordType P0(Mesh::ConstFacePointer fp, int i = 0);
    Mesh::CoordType P1(Mesh::ConstFacePointer fp, int i = 0);
    Mesh::CoordType P2(Mesh::ConstFacePointer fp, int i = 0);

};


class SymmetricDirichlet : public Energy {

    SimpleTempData<Mesh::FaceContainer, Point4d> data;

public:

    SymmetricDirichlet(Mesh& mesh, Geometry geometryMode);

    double E(const Mesh::FaceType& f);

    Eigen::MatrixXd Grad();
};


#endif // ENERGY_H

