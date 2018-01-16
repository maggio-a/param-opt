#ifndef ITERATIVE_H
#define ITERATIVE_H

#include "mesh.h"
#include "energy.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <memory>

class DescentMethod {

protected:

    Mesh& m;
    std::shared_ptr<Energy> energy;

public:

    DescentMethod(std::shared_ptr<Energy> e) : m{e->m}, energy{e} {}
    virtual ~DescentMethod() {}

    virtual Eigen::MatrixXd ComputeDescentDirection() = 0;

    double SearchStrongWolfe(const Eigen::MatrixXd& uv, const Eigen::MatrixXd& grad, const Eigen::MatrixXd& dir);
    double Search(const Eigen::MatrixXd& uv, const Eigen::MatrixXd& grad, const Eigen::MatrixXd& dir);
    virtual double Iterate(double& gradientNorm, double& objValDiff);

    Eigen::MatrixXd X();
    void SetX(const Eigen::MatrixXd& x);
};


class GradientDescent : public DescentMethod {

public:

    GradientDescent(std::shared_ptr<Energy> energy) : DescentMethod(energy) {}

    virtual Eigen::MatrixXd ComputeDescentDirection();
};


class LBFGS : public DescentMethod {

    const std::size_t m;
    std::deque<Eigen::MatrixXd> sVec;
    std::deque<Eigen::MatrixXd> yVec;
    std::deque<double> rho;

public:

    LBFGS(std::shared_ptr<Energy> energy, std::size_t memory) : DescentMethod(energy), m{memory} {}

    Eigen::MatrixXd ComputeDescentDirection();
    virtual double Iterate(double& gradientNorm, double& objValDiff);
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

public:

    SLIM(std::shared_ptr<SymmetricDirichlet> sd);

    virtual Eigen::MatrixXd ComputeDescentDirection();

private:

    void UpdateJRW();
    void BuildA(Eigen::SparseMatrix<double>& A);
    void BuildRhs(const Eigen::SparseMatrix<double>& At, Eigen::VectorXd &rhs);
    Eigen::MatrixXd MinimizeProxyEnergy();
};

#endif // ITERATIVE_H

