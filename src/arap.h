#ifndef ARAP_H
#define ARAP_H

#include "mesh.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

class ARAP {

public:

    struct Cot {
        double v[3];
    };

private:

    Mesh& m;

    std::vector<int> fixed_i;
    std::vector<vcg::Point2d> fixed_pos;

    int max_iter;

    double CurrentEnergy();
    void ComputeSystemMatrix(Mesh& m, const std::vector<Cot>& cotan, Eigen::SparseMatrix<double>& L);
    void ComputeRHS(Mesh& m, const std::vector<Eigen::Matrix2d>& rotations, const std::vector<Cot>& cotan, Eigen::VectorXd& bu, Eigen::VectorXd& bv);

public:

    ARAP(Mesh& mesh);

    void FixVertex(Mesh::ConstVertexPointer vp, const vcg::Point2d& pos);
    void FixBoundaryVertices();
    int FixSelectedVertices();
    int FixRandomEdgeWithinTolerance(double tol);
    void SetMaxIterations(int n);

    void Solve();
};



#endif // ARAP_H
