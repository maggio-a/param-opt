#ifndef ENERGY_H
#define ENERGY_H

#include "mesh.h"
#include "uv.h"
#include "mesh_attribute.h"

#include <Eigen/Core>

/*
 *
 * A note on the management of scaffold triangles. The issue with the scaffold
 * is that for a non-trivial number of iterations on meshes that overlap, the
 * scaffold triangle will degenerate to slivers. This seems to cause issues with
 * the computation of the descent direction using SLIM, so the idea is to set a
 * threshold value on the quality of the scaffold triangles, below which the
 * energy contribution of the scaffold triangles is not dampened. This should
 * also avoid generating parameterizations where the boundary of the chart gets
 * arbitrarily close to itself without intersecting.
 *
 * In order to do so, an index of the quality of scaffold triangles is cached in
 * the energy object, which is used to compute the actual energy contributions.
 * Let q be the quality threshold. For a given scaffold face s_f, if Q[s_f] < q
 * the energy contribution of s_f behaves exactly as any other face, otherwise
 * the contribution is normalized according to the SCAF paper, changing the
 * total number of scaffold faces in the energy term with the number of faces
 * that are above the quality threshold.
 *
 * It turns out that simply treating the near-degenerate s_f as a mesh face is
 * not enough to prevent slivers, because the rest pose of the face is exactly
 * the starting configuration of the face, and therefore the optimizer induces
 * small distortions that still shrink those faces, eventually yielding
 * extremely close boundaries that Triangle triangulates with slivers. The
 * solution implemented here heavily penalizes deformations applied to faces
 * that are near-degenerate, essentially locking them. This approach has the
 * obvious defect of preventing further optimization of the mesh faces nearby.
 *
 */

enum EnergyType {
    SymmetricDirichlet
};

class Energy {

    friend class DescentMethod;
    friend class SLIM;
    friend class CompMaj;

protected:

    Mesh& m;
    double surfaceArea;
    double holeFillingArea;
    Mesh::PerFaceAttributeHandle<CoordStorage> targetShape;

    int numMeshFaces;
    int numHoleFillingFaces;
    int numScaffoldFaces;

    double scaffoldRegularizationTerm;

    virtual double E(const Mesh::FaceType& f) = 0;
    virtual void Grad(int faceIndex, Eigen::Vector2d& g0, Eigen::Vector2d& g1, Eigen::Vector2d& g2) = 0;

public:

    Energy(Mesh& mesh);
    virtual ~Energy();

    double E();
    double E_IgnoreMarkedFaces();

    Eigen::VectorXd EnergyPerFace();

    virtual double NormalizedMinValue() = 0;

    Eigen::MatrixXd Grad();

    double GetScaffoldWeight() const;
    double GetScaffoldFaceArea() const;

    /* Utility function to update cached data. This must be called whenever the
     * underlying mesh changes. We need it to update the cached data of the
     * hole-filling faces when they are remeshed.
     * In order to update the Energy obj cache with the scaffold weight, I need
     * to be able to compute the energy value in the current configuration, but
     * to do so the actual energy object needs to have its cacheed data already
     * updated. This means that subclasses are required to call Energy::UpdateCache()
     * only _after_ their cache has already been updated. */
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

inline Mesh::CoordType Energy::P(Mesh::ConstFacePointer fp, int i)
{
   assert(i >= 0 && i <= 2);
   return targetShape[fp].P[i];
}

inline Mesh::CoordType Energy::P0(Mesh::ConstFacePointer fp, int i) { return P(fp, i); }
inline Mesh::CoordType Energy::P1(Mesh::ConstFacePointer fp, int i) { return P(fp, (i+1)%3); }
inline Mesh::CoordType Energy::P2(Mesh::ConstFacePointer fp, int i) { return P(fp, (i+2)%3); }

inline double Energy::GetScaffoldWeight() const
{
    assert(scaffoldRegularizationTerm > 0);
    return scaffoldRegularizationTerm;
}

inline double Energy::GetScaffoldFaceArea() const
{
    return 1.0;
}


class SymmetricDirichletEnergy : public Energy {

    friend class SLIM;
    friend class CompMaj;

    /* Precomputed cotangents and area of each target shape */
    SimpleTempData<Mesh::FaceContainer, Point4d> data;

public:

    SymmetricDirichletEnergy(Mesh& mesh);

    double NormalizedMinValue();
    void UpdateCache();

protected:

    double E(const Mesh::FaceType& f) override;
    void Grad(int faceIndex, Eigen::Vector2d& g0, Eigen::Vector2d& g1, Eigen::Vector2d& g2) override;
};


#endif // ENERGY_H

