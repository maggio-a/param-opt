#ifndef ITERATIVE_H
#define ITERATIVE_H

#include "mesh.h"
#include "energy.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <memory>

enum DescentType {
    Gradient, LimitedMemoryBFGS, ScalableLocallyInjectiveMappings
};

class DescentMethod {

protected:

    Mesh& m;
    std::shared_ptr<Energy> energy;

public:

    DescentMethod(std::shared_ptr<Energy> e);
    virtual ~DescentMethod();

    // Interface

    virtual Eigen::MatrixXd ComputeDescentDirection() = 0;
    virtual double Iterate(double& gradientNorm, double& objValDiff);

    /* Utility function to update cached data. This must be called whenever the
     * underlying mesh changes. We need it to update the cached data of the
     * hole-filling faces when they are remeshed. */
    virtual void UpdateCache();


    // Member functions

    /* Line search that satisfies the strong Wolfe conditions.
     * Reference: Nocedal and Wright (2006), "Numerical Optimization" */
    double SearchStrongWolfe(const Eigen::MatrixXd& uv, const Eigen::MatrixXd& grad, const Eigen::MatrixXd& dir);

    /* Simple bisection based line search (Armijo) */
    double Search(const Eigen::MatrixXd& uv, const Eigen::MatrixXd& grad, const Eigen::MatrixXd& dir);

    /* Returns the current solution as a VN-by-2 matrix of doubles */
    Eigen::MatrixXd X();

    /* Sets VN-by-2 matrix of doubles as the vertex texture coordinates of the
     * underlying mesh object */
    void SetX(const Eigen::MatrixXd& x);
};


class GradientDescent : public DescentMethod {

public:

    GradientDescent(std::shared_ptr<Energy> energy);

    virtual Eigen::MatrixXd ComputeDescentDirection();
};


class LBFGS : public DescentMethod {

    const std::size_t m;
    std::deque<Eigen::MatrixXd> sVec;
    std::deque<Eigen::MatrixXd> yVec;
    std::deque<double> rho;

public:

    LBFGS(std::shared_ptr<Energy> energy, std::size_t memory);

    Eigen::MatrixXd ComputeDescentDirection();
    virtual double Iterate(double& gradientNorm, double& objValDiff);
    virtual void UpdateCache();
};


/*
 * Implementation of Rabinovich et al. (2017) 'Scalable Locally Injective Mappings'
 */
class SLIM : public DescentMethod {

    SimpleTempData<Mesh::FaceContainer, Eigen::Matrix2d> fm_inv; // 2by2 map from mesh face to canonical triangle
    SimpleTempData<Mesh::FaceContainer, Eigen::Matrix2d> J; // per face jacobians
    SimpleTempData<Mesh::FaceContainer, Eigen::Matrix2d> R; // per face closest rotations
    SimpleTempData<Mesh::FaceContainer, Eigen::Matrix2d> W; // per face weight coefficients

    Eigen::SparseMatrix<double> D1; // Dx in paper appendix
    Eigen::SparseMatrix<double> D2; // Dy in paper appendix
    Eigen::VectorXd diagAreaVector; // Vector of face areas repeated 4 times (used as diagonal matrix)
    const double lambda; // the 'proximal penalty' term in the paper

    // store the solver for the global step in order to reuse the factorization pattern
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    bool firstSolve;

public:

    SLIM(std::shared_ptr<SymmetricDirichletEnergy> sd);

    virtual Eigen::MatrixXd ComputeDescentDirection();
    virtual void UpdateCache();

private:

    void PrecomputeD12(const Mesh& m, Eigen::SparseMatrix<double>& D1, Eigen::SparseMatrix<double>& D2);
    void UpdateJRW();
    void BuildA(Eigen::SparseMatrix<double>& A);
    void BuildRhs(const Eigen::SparseMatrix<double>& At, Eigen::VectorXd &rhs);
    void MinimizeProxyEnergy(Eigen::MatrixXd &p_k);
};

#endif // ITERATIVE_H

