#ifndef ENERGY_H
#define ENERGY_H

#include "mesh.h"
#include "uv.h"
#include "mesh_attribute.h"

#include <Eigen/Core>


class Energy {

    friend class DescentMethod;
    friend class SLIM;

protected:

    Mesh& m;
    double surfaceArea;
    double holeFillingArea;
    Mesh::PerFaceAttributeHandle<CoordStorage> targetShape;

public:

    Energy(Mesh& mesh);
    virtual ~Energy();

    double E();
    double E_IgnoreMarkedFaces(bool normalized = false);
    virtual double E(const Mesh::FaceType& f, bool normalized = false) = 0;
    virtual Eigen::MatrixXd Grad() = 0;

    virtual void UpdateCache();

    void MapToFaceQuality(bool normalized);

    double FaceArea(Mesh::ConstFacePointer fp);
    double SurfaceArea();
    double ParameterArea();

    /* Hack to speed up the convergence of iterative methods
     * Scales the parameterization so that its area matches the total area of the model (either the actual 3D area or
     * the area of the original parameterized faces, according to the geometry mode).
     * This can be used by iterative methods on energies that penalize area mismatches, in order to scale the initial solution
     * and reduce the scaling process that happens during the iterations (thereby lowering the number of iterations required
     * to converge) */
    void CorrectScale();

    Mesh::CoordType P(Mesh::ConstFacePointer fp, int i);
    Mesh::CoordType P0(Mesh::ConstFacePointer fp, int i = 0);
    Mesh::CoordType P1(Mesh::ConstFacePointer fp, int i = 0);
    Mesh::CoordType P2(Mesh::ConstFacePointer fp, int i = 0);

};

inline Mesh::CoordType Energy::P0(Mesh::ConstFacePointer fp, int i) { return P(fp, i); }
inline Mesh::CoordType Energy::P1(Mesh::ConstFacePointer fp, int i) { return P(fp, (i+1)%3); }
inline Mesh::CoordType Energy::P2(Mesh::ConstFacePointer fp, int i) { return P(fp, (i+2)%3); }

class SymmetricDirichlet : public Energy {

    SimpleTempData<Mesh::FaceContainer, Point4d> data;

public:

    SymmetricDirichlet(Mesh& mesh);

    double E(const Mesh::FaceType& f, bool normalized = false);
    Eigen::MatrixXd Grad();
    void UpdateCache();
};


#endif // ENERGY_H

