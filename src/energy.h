#ifndef ENERGY_H
#define ENERGY_H

#include "mesh.h"
#include "uv.h"
#include "mesh_attribute.h"

#include <Eigen/Core>

enum EnergyType {
    SymmetricDirichlet
};

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

    virtual double NormalizedMinValue() = 0;

    Eigen::MatrixXd Grad();
    virtual void Grad(int faceIndex, Eigen::Vector2d& g0, Eigen::Vector2d& g1, Eigen::Vector2d& g2) = 0;

    /* Utility function to update cached data. This must be called whenever the
     * underlying mesh changes. We need it to update the cached data of the
     * hole-filling faces when they are remeshed. */
    virtual void UpdateCache();

    void MapToFaceQuality(bool normalized);

    double FaceArea(Mesh::ConstFacePointer fp);
    double SurfaceArea();
    double SurfaceAreaNotHoleFilling();
    double ParameterArea();

    /* Hack to speed up the convergence of iterative methods.
     * Scales the parameterization so that its area matches the total area of
     * the model (either the actual 3D area or the area of the original
     * parameterized faces, according to the geometry mode). This can be used by
     * iterative methods on energies that penalize area mismatches, in order to
     * scale the initial solution and reduce the scaling process that happens
     * during the iterations (thereby lowering the number of iterations required
     * to converge) */
    void CorrectScale();

    /* Accessor functions to retrieve the positions of the target shape for each
     * face */
    Mesh::CoordType P(Mesh::ConstFacePointer fp, int i);
    Mesh::CoordType P0(Mesh::ConstFacePointer fp, int i = 0);
    Mesh::CoordType P1(Mesh::ConstFacePointer fp, int i = 0);
    Mesh::CoordType P2(Mesh::ConstFacePointer fp, int i = 0);

};

inline Mesh::CoordType Energy::P0(Mesh::ConstFacePointer fp, int i) { return P(fp, i); }
inline Mesh::CoordType Energy::P1(Mesh::ConstFacePointer fp, int i) { return P(fp, (i+1)%3); }
inline Mesh::CoordType Energy::P2(Mesh::ConstFacePointer fp, int i) { return P(fp, (i+2)%3); }

class SymmetricDirichletEnergy : public Energy {

    /* Precomputed cotangents and area of each target shape */
    SimpleTempData<Mesh::FaceContainer, Point4d> data;

public:

    SymmetricDirichletEnergy(Mesh& mesh);

    double E(const Mesh::FaceType& f, bool normalized = false);
    double NormalizedMinValue();
    void Grad(int faceIndex, Eigen::Vector2d& g0, Eigen::Vector2d& g1, Eigen::Vector2d& g2);
    void UpdateCache();
};


#endif // ENERGY_H

