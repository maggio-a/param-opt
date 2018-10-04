#ifndef ITERATIVE_H
#define ITERATIVE_H

#include "mesh.h"
#include "energy.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <memory>

using Matrix66d = Eigen::Matrix<double, 6, 6>;
using Vector6d  = Eigen::Matrix<double, 6, 1>;

enum DescentType {
    Gradient, LimitedMemoryBFGS, ScalableLocallyInjectiveMappings, CompositeMajorization
};

class DescentMethod {

protected:

    Mesh& m;
    std::shared_ptr<Energy> energy;

public:

    DescentMethod(std::shared_ptr<Energy> e);
    virtual ~DescentMethod();

    // Interface

    virtual bool ComputeDescentDirection(Eigen::MatrixXd& dir) = 0;
    virtual bool Iterate(double& gradientNorm, double& objValDiff, double& energyAfterIter);

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

    virtual bool ComputeDescentDirection(Eigen::MatrixXd& dir);
};


class LBFGS : public DescentMethod {

    const std::size_t m;
    std::deque<Eigen::MatrixXd> sVec;
    std::deque<Eigen::MatrixXd> yVec;
    std::deque<double> rho;

public:

    LBFGS(std::shared_ptr<Energy> energy, std::size_t memory);

    virtual bool ComputeDescentDirection(Eigen::MatrixXd& dir);
    virtual bool Iterate(double& gradientNorm, double& objValDiff, double& energyAfterIter);
    virtual void UpdateCache();
};


/*
 * Rabinovich et al. (2017) 'Scalable Locally Injective Mappings'
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

    std::shared_ptr<SymmetricDirichletEnergy> sd_energy;

public:

    SLIM(std::shared_ptr<SymmetricDirichletEnergy> sd);

    virtual bool ComputeDescentDirection(Eigen::MatrixXd& dir);
    virtual void UpdateCache();

private:

    void PrecomputeD12(const Mesh& m, Eigen::SparseMatrix<double>& D1, Eigen::SparseMatrix<double>& D2);
    void UpdateJRW();
    void BuildA(Eigen::SparseMatrix<double>& A);
    void BuildRhs(const Eigen::SparseMatrix<double>& At, Eigen::VectorXd &rhs);
    bool MinimizeProxyEnergy(Eigen::MatrixXd &p_k);
    void BuildSystem(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &rhs);
};


/*
 * Shtengel et al. (2017) 'Geometric optimization via Composite Majorization'
 */
class CompMaj : public DescentMethod {

public:

    CompMaj(std::shared_ptr<SymmetricDirichletEnergy> sd);

    virtual bool ComputeDescentDirection(Eigen::MatrixXd& dir);
    virtual void UpdateCache();

private:

    std::shared_ptr<SymmetricDirichletEnergy> sd_energy;

    int vn;
    int fn;

    Eigen::VectorXd a;
    Eigen::VectorXd b;
    Eigen::VectorXd c;
    Eigen::VectorXd d;

    Eigen::VectorXd detJuv;

    Eigen::MatrixX2d s;
    Eigen::MatrixX4d v;
    Eigen::MatrixX4d u;
    Eigen::MatrixXd Dsd[2];

    Eigen::Matrix3Xd D1d;
    Eigen::Matrix3Xd D2d;

    Eigen::MatrixXd a1d;
    Eigen::MatrixXd a2d;
    Eigen::MatrixXd b1d;
    Eigen::MatrixXd b2d;

    std::vector<Matrix66d> Hi; // vector of hessian matrices per face

    std::vector<int> II;
    std::vector<int> JJ;
    std::vector<double> SS;

    Eigen::SparseMatrix<double> A;

    void ComputeSurfaceGradientPerFace(Eigen::MatrixX3d& D1, Eigen::MatrixX3d& D2);

    void UpdateJ();
    void ComputeHessian();
    void UpdateSSVDFunction();
    void ComputeDenseSSVDDerivatives();

    Matrix66d ComputeFaceConeHessian(const Vector6d& A1, const Vector6d& A2, double a1x, double a2x);

    Matrix66d ComputeConvexConcaveFaceHessian(const Vector6d& a1, const Vector6d& a2, const Vector6d& b1, const Vector6d& b2,
            double alpha0, double alpha1, double beta0, double beta1,
            const Vector6d& dS0i, const Vector6d& dS1i,
            double gradS0, double gradS1, double HS0, double HS1);

};

#endif // ITERATIVE_H
