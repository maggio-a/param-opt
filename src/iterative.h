#ifndef ITERATIVE_H
#define ITERATIVE_H

#include "mesh.h"
#include "energy.h"

#include <Eigen/Core>

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

    Eigen::MatrixXd ComputeDescentDirection();

};


class LBFGS : public DescentMethod {

    const std::size_t m;
    std::deque<Eigen::MatrixXd> sVec;
    std::deque<Eigen::MatrixXd> yVec;
    std::deque<double> rho;

public:

    LBFGS(std::shared_ptr<Energy> energy, std::size_t memory) : DescentMethod(energy), m{memory} {}

    Eigen::MatrixXd ComputeDescentDirection();
    double Iterate(double& gradientNorm, double& objValDiff);
};

#endif // ITERATIVE_H

