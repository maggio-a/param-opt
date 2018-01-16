#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "iterative.h"
#include "mesh.h"
#include "energy.h"
#include "math_utils.h"
#include "metric.h"

using Eigen::MatrixXd;
using Eigen::Vector2d;

// Quadratic equation solver
static double QuadraticRoots(double a, double b, double c, double *x1, double *x2)
{
    double delta = b*b - 4.0f*a*c;
    if (delta >= 0) {
        *x1 = (-b + std::sqrt(delta)) / (2*a);
        *x2 = (-b - std::sqrt(delta)) / (2*a);
    }
    return delta;
}

// this functions takes the 2d base coordinates and the offsets for each vertex, and returns the smallest positive
// step that prevents a triangle flip, or numeric_limits::max() if the step is unbounded
// Negli appunti pij = (u_ij, v_ij), uguale per pki
static double ComputeStepSizeNoFlip(const Vector2d& pi, const Vector2d& pj, const Vector2d& pk, const Vector2d& di, const Vector2d& dj, const Vector2d& dk)
{
    Vector2d dji = dj - di; Vector2d dki = dk - di;
    Vector2d pji = pj - pi; Vector2d pki = pk - pi;

    double a = dji[0]*dki[1] - dji[1]*dki[0];
    double b = pji[0]*dki[1] + pki[1]*dji[0] - pji[1]*dki[0] - pki[0]*dji[1];
    double c = pji[0]*pki[1] - pji[1]*pki[0];

    double x1, x2;
    double delta = QuadraticRoots(a, b, c, &x1, &x2);
    if (delta < 0) return std::numeric_limits<double>::max();
    //assert(x1 != 0 && x2 != 0);
    if (x2 < x1) std::swap(x1, x2);
    if (x1 > 0) return x1;
    else if (x2 > 0) return x2;
    else return std::numeric_limits<double>::max();
}


// Descent Method implementation
// =============================

double DescentMethod::SearchStrongWolfe(const MatrixXd& uv, const MatrixXd& grad, const MatrixXd& dir)
{
    double c1 = 0.2;
    double c2 = 0.9;
    double phi_zero = energy->E();
    double phiprime_zero = dir.cwiseProduct(grad).sum();

    double maxStepWithoutInversion = std::numeric_limits<double>::max();
    for (auto& f : m.face) {
        //assert(((f.V(1)->T().P() - f.V(0)->T().P()) ^ (f.V(2)->T().P() - f.V(0)->T().P())) > 0);
        double localStep = ComputeStepSizeNoFlip(uv.row(Index(m, f.V(0))), uv.row(Index(m, f.V(1))), uv.row(Index(m, f.V(2))),
                                                 dir.row(Index(m, f.V(0))), dir.row(Index(m, f.V(1))), dir.row(Index(m, f.V(2))));
        if (localStep < maxStepWithoutInversion) maxStepWithoutInversion = localStep;
        /*float t = std::min(localStep, 10000.0); //std::min(0.8*localStep, 100000.0);
        Vector2d a = uv.row(Index(m, f.V(0))) + t * dir.row(Index(m, f.V(0)));
        Vector2d b = uv.row(Index(m, f.V(1))) + t * dir.row(Index(m, f.V(1)));
        Vector2d c = uv.row(Index(m, f.V(2))) + t * dir.row(Index(m, f.V(2)));
        Point2d d1{b[0] - a[0], b[1] - a[1]};
        Point2d d2{c[0] - a[0], c[1] - a[1]};
        double aaa = d1 ^ d2;
        if (aaa <= 0) {
            std::cout << t << "   " << aaa << std::endl;
            Vector2d u0 = uv.row(Index(m, f.V(0)));
            Vector2d u1 = uv.row(Index(m, f.V(1)));
            Vector2d u2 = uv.row(Index(m, f.V(2)));

            Vector2d d0 = dir.row(Index(m, f.V(0)));
            Vector2d d1 = dir.row(Index(m, f.V(1)));
            Vector2d d2 = dir.row(Index(m, f.V(2)));

            std::cout << u0[0] << " " << u0[1] << " dir =  " << d0[0] << " " << d0[1] << std::endl;
            std::cout << u1[0] << " " << u1[1] << " dir =  " << d1[0] << " " << d1[1] << std::endl;
            std::cout << u2[0] << " " << u2[1] << " dir =  " << d2[0] << " " << d2[1] << std::endl;

        }
        assert ((d1^d2) >= 0);*/
    }

    double alpha_max = std::min(1.0 / 0.8, 0.99 * maxStepWithoutInversion);
    double alpha_prev = 0;
    double alpha_curr = 0.8 * alpha_max;

    double phi_prev = phi_zero;

    double alpha_lo = std::numeric_limits<double>::lowest();
    double alpha_hi = std::numeric_limits<double>::lowest();
    double phi_alpha_lo = 0;

    int i = 1;
    while (alpha_curr < alpha_max) {
        SetX(uv + alpha_curr * dir);
        double phi_curr = energy->E();
        if ((phi_curr > phi_zero + c1 * alpha_curr * phiprime_zero) || (phi_prev >= phi_curr && i > 1)) {
            alpha_lo = alpha_prev;
            phi_alpha_lo = phi_prev;
            alpha_hi = alpha_curr;
            break;
        }
        double phiprime_curr = dir.cwiseProduct(energy->Grad()).sum();
        if (std::abs(phiprime_curr) <= c2 * phiprime_zero) {
            assert(alpha_curr < maxStepWithoutInversion);
            return alpha_curr;
        }
        if (phiprime_curr >= 0) {
            alpha_lo = alpha_curr;
            phi_alpha_lo = phi_curr;
            alpha_hi = alpha_prev;
            break;
        }
        // setup next iteration
        //double alpha_next = (alpha_curr + alpha_max) / 2.0;
        alpha_prev = alpha_curr;
        phi_prev = phi_curr;
        alpha_curr = (alpha_curr + alpha_max) / 2.0;
        ++i;
    }

    if (alpha_curr >= alpha_max) {
        return alpha_max;
    }

    int k = 0;
    constexpr int MAX_ITERATIONS = 10;
    while (k < MAX_ITERATIONS) {
        alpha_curr = (alpha_lo + alpha_hi) / 2.0;
        SetX(uv + alpha_curr * dir);
        double phi_curr = energy->E();
        if ((phi_curr > phi_zero + c1 * alpha_curr * phiprime_zero) || (phi_curr >= phi_alpha_lo)) {
            alpha_hi = alpha_curr;
        }
        else {
            double phiprime_curr = dir.cwiseProduct(energy->Grad()).sum();
            if (std::abs(phiprime_curr) <= -c2 * phiprime_zero) {
                assert(alpha_curr < maxStepWithoutInversion);
                return alpha_curr;
            }
            if (phiprime_curr * (alpha_hi - alpha_lo) >= 0) {
                alpha_hi = alpha_lo;
            }
            alpha_lo = alpha_curr;
            phi_alpha_lo = phi_curr;
        }
        k++;
    }

    assert(alpha_curr < maxStepWithoutInversion);
    return alpha_curr;
}

double DescentMethod::Search(const MatrixXd& uv, const MatrixXd& grad, const MatrixXd& dir)
{
    double E_curr = energy->E();

    double c = 0.2;
    double descentCoeff = c * dir.cwiseProduct(grad).sum();

    assert(descentCoeff < 0); // otherwise dir is not a descent direction

    double t = std::numeric_limits<double>::max();
    for (auto& f : m.face) {
        double tFace = ComputeStepSizeNoFlip(uv.row(Index(m, f.V(0))),
                                             uv.row(Index(m, f.V(1))),
                                             uv.row(Index(m, f.V(2))),
                                             dir.row(Index(m, f.V(0))),
                                             dir.row(Index(m, f.V(1))),
                                             dir.row(Index(m, f.V(2))));
        if (tFace < t) t = tFace;
    }
    t = std::min(1.0, t*0.8);

    SetX(uv + (t * dir));
    double E_t = energy->E();
    int numIter = 0;
    double gamma = 0.8;
    while (E_t > E_curr + t * descentCoeff) {
        t *= gamma;
        SetX(uv + (t * dir));
        E_t = energy->E();
        numIter++;
    }

    std::cout << "t = " << t << std::endl << "descentCoeff = " << descentCoeff << std::endl;
    std::cout << "Descent took " << numIter << " iterations to decrease function" << std::endl;

    return E_t;
}

double DescentMethod::Iterate(double& gradientNorm, double& objValDiff)
{
    double energyPrev = energy->E();

    MatrixXd uv = X();
    MatrixXd grad = energy->Grad();
    MatrixXd dir = ComputeDescentDirection();
    Search(uv, grad, dir);

    double energyCurr = energy->E();
    MatrixXd newGrad = energy->Grad();
    gradientNorm = std::sqrt(newGrad.cwiseProduct(newGrad).sum());
    objValDiff = energyPrev - energyCurr;

    return energyCurr;
}

MatrixXd DescentMethod::X()
{
    MatrixXd x{m.VN(), 2};

    for (auto& v : m.vert) {
        x.row(Index(m, v)) = Eigen::Vector2d{v.T().U(), v.T().V()};
    }

    return x;
}

void DescentMethod::SetX(const MatrixXd& x)
{
    assert(x.rows() == m.VN() && x.cols() == 2);
    for (auto& v : m.vert) {
        auto uv = x.row(Index(m, v));
        v.T().U() = uv(0);
        v.T().V() = uv(1);
    }
}

// Gradient Descent
// ================

MatrixXd GradientDescent::ComputeDescentDirection()
{
    return - energy->Grad();
}

// L-BFGS
// ======

MatrixXd LBFGS::ComputeDescentDirection()
{
    assert(sVec.size() == yVec.size());
    int memsz = sVec.size();
    MatrixXd q = energy->Grad();
    if (memsz > 0) {
        std::vector<double> alpha;
        alpha.reserve(memsz);
        for (int i = 0; i < memsz; ++i) {
            alpha.push_back(rho[i] * sVec[i].cwiseProduct(q).sum());
            q = q - alpha[i] * yVec[i];
        }
        //q *= (sVec[0].cwiseProduct(yVec[0]).sum()) / (yVec[0].cwiseProduct(yVec[0]).sum());
        q *= (1.0 / rho[0]) / (yVec[0].cwiseProduct(yVec[0]).sum());
        for (int i = memsz - 1; i >= 0; --i) {
            double b_i = rho[i] * yVec[i].cwiseProduct(q).sum();
            q = q + (alpha[i] - b_i) * sVec[i];
        }
    }
    return -q;
}

double LBFGS::Iterate(double& gradientNorm, double& objValDiff)
{
    //energy->CorrectScale();
    double energyPrev = energy->E();
    MatrixXd uv = X();
    MatrixXd grad = energy->Grad();
    MatrixXd dir = ComputeDescentDirection();
    //MatrixXd sdir = ComputeScaleDirection();
    //MatrixXd combined = dir + sdir;
    //if (combined.cwiseProduct(grad).sum() < 0) {
        //std::cout << "Using combined dir" << std::endl;
    //    dir = combined;
    //}
    SearchStrongWolfe(uv, grad, dir);

    double energyCurr = energy->E();
    MatrixXd newGrad = energy->Grad();

    // manage lbfgs state
    sVec.push_front(X() - uv);
    yVec.push_front(newGrad - grad);
    rho.push_front(1.0 / yVec[0].cwiseProduct(sVec[0]).sum());
    if (sVec.size() > m) sVec.pop_back();
    if (yVec.size() > m) yVec.pop_back();
    if (rho.size() > m) rho.pop_back();

    gradientNorm = std::sqrt(newGrad.cwiseProduct(newGrad).sum());
    objValDiff = energyPrev - energyCurr;

    return energyCurr;
}

// Scalable locally injective mappings
// ===================================


/*
 * CREDITS
 *
 * The code for the least squares solver is based on the reference implemenation from libigl available at
 *
 * https://github.com/libigl/libigl
 *
 * and in particular the file slim.cpp
 *
 * */

static void PrecomputeD12(const Mesh& m, Eigen::SparseMatrix<double>& D1, Eigen::SparseMatrix<double>& D2);

SLIM::SLIM(std::shared_ptr<SymmetricDirichlet> sd)
    : DescentMethod(sd),
      fm_inv{sd->m.face},
      J{sd->m.face},
      R{sd->m.face},
      W{sd->m.face},
      D1{},
      D2{},
      diagAreaVector{},
      lambda{0.0001}
{
    for (auto& f : m.face) {
        /*
         * find 2d coords for the vertices of the faces relative to the frame which is coplanar with
         * the face and has the origin in P0. x0 is (0,0), x1 is (lenght(p1-p0), 0) and x2 is (1,0) rotated
         * by the angle betwen (p1-p0) and (p2-p0)
         * */
        Point3d p10 = energy->P1(&f) - energy->P0(&f);
        Point3d p20 = energy->P2(&f) - energy->P0(&f);
        double theta = VecAngle(p10, p20);
        Eigen::Vector2d x1{p10.Norm(), 0};
        double st = std::sin(theta);
        double ct = std::cos(theta);
        Eigen::Matrix2d rt;
        rt << ct, -st, st, ct;
        Eigen::Vector2d x2 = rt * Eigen::Vector2d{p20.Norm(), 0}; // rotation after scaling

        /*
         * at this point the transformation from the canonical triangle (0,0) (1,0), (0,1) to the
         * mesh face is encoded in the matrix [x1 | x2], and we want the inverse
         * inverse = (1/det) * adjugate
         * */
        double a = x1(0), b = x2(0), c = x1(1), d = x2(1);
        fm_inv[f] << d, -b, -c, a;
        fm_inv[f] *= 1.0 / (b*d - a*c);
    }

    PrecomputeD12(m, D1, D2);

    diagAreaVector.resize(4 * m.FN());
    for (auto const& f : m.face) {
        double area = 0;
        switch (energy->mode) {
        case Model: area = DistortionMetric::Area3D(m, f, ParameterizationGeometry::Model); break;
        case Texture: area = DistortionMetric::Area3D(m, f, ParameterizationGeometry::Texture); break;
        }

        int j = tri::Index(m, f);
        diagAreaVector(0 * m.FN() + j) = area;
        diagAreaVector(1 * m.FN() + j) = area;
        diagAreaVector(2 * m.FN() + j) = area;
        diagAreaVector(3 * m.FN() + j) = area;
    }
}

static void PrecomputeD12(const Mesh& m, Eigen::SparseMatrix<double>& D1, Eigen::SparseMatrix<double>& D2)
{
    Eigen::Matrix<double, Eigen::Dynamic, 3> eperp21(m.FN(), 3), eperp13(m.FN(), 3);
    Eigen::Matrix<double, Eigen::Dynamic, 3> F1(m.FN(), 3), F2(m.FN(), 3); // basis vectors for the tangent plane at each face

    //for (int i=0;i<F.rows();++i)
    for (auto const& f : m.face) {
        // #F x 3 matrices of triangle edge vectors, named after opposite vertices
        Eigen::Matrix<double, 1, 3> v32; (f.cP(2) - f.cP(1)).ToEigenVector(v32); // V.row(i3) - V.row(i2);
        Eigen::Matrix<double, 1, 3> v13; (f.cP(0) - f.cP(2)).ToEigenVector(v13); // V.row(i1) - V.row(i3);
        Eigen::Matrix<double, 1, 3> v21; (f.cP(1) - f.cP(0)).ToEigenVector(v21); // V.row(i2) - V.row(i1);
        Eigen::Matrix<double, 1, 3> n = v32.cross(v13);
        // area of parallelogram is twice area of triangle
        // area of parallelogram is || v1 x v2 ||
        // This does correct l2 norm of rows, so that it contains #F list of twice
        // triangle areas
        double dblA = std::sqrt(n.dot(n));
        Eigen::Matrix<double, 1, 3> u = n / dblA;

        int i = tri::Index(m, f); // face index

        F1.row(i) = v21.normalized();
        F2.row(i) = u.cross(v21).normalized();

        // rotate each vector 90 degrees around normal
        double norm21 = std::sqrt(v21.dot(v21));
        double norm13 = std::sqrt(v13.dot(v13));
        //eperp21.row(i) = u.cross(v21);
        //eperp21.row(i) = eperp21.row(i) / std::sqrt(eperp21.row(i).dot(eperp21.row(i)));
        eperp21.row(i) = F2.row(i);
        eperp21.row(i) *= norm21 / 0.5*dblA;
        eperp13.row(i) = u.cross(v13);
        eperp13.row(i) = eperp13.row(i) / std::sqrt(eperp13.row(i).dot(eperp13.row(i)));
        eperp13.row(i) *= norm13 / dblA;
    }

    std::vector<int> rs;
    rs.reserve(m.FN()*4*3);
    std::vector<int> cs;
    cs.reserve(m.FN()*4*3);
    std::vector<double> vs;
    vs.reserve(m.FN()*4*3);

    // row indices
    for(int r = 0; r < 3; r++) {
        for(int j = 0; j < 4; j++) {
            for(int i = r*m.FN(); i < (r+1)*m.FN(); i++) rs.push_back(i);
        }
    }

    // column indices
    for (int r = 0; r < 3; r++) {
        for (int i = 0; i < m.FN(); i++) cs.push_back(tri::Index(m, m.face[i].cV(1)));//F(i,1));
        for (int i = 0; i < m.FN(); i++) cs.push_back(tri::Index(m, m.face[i].cV(0)));//F(i,0));
        for (int i = 0; i < m.FN(); i++) cs.push_back(tri::Index(m, m.face[i].cV(2)));//F(i,2));
        for (int i = 0; i < m.FN(); i++) cs.push_back(tri::Index(m, m.face[i].cV(0)));//F(i,0));
    }

    // values
    for (int i = 0; i < m.FN(); i++) vs.push_back(eperp13(i,0));
    for (int i = 0; i < m.FN(); i++) vs.push_back(-eperp13(i,0));
    for (int i = 0; i < m.FN(); i++) vs.push_back(eperp21(i,0));
    for (int i = 0; i < m.FN(); i++) vs.push_back(-eperp21(i,0));
    for (int i = 0; i < m.FN(); i++) vs.push_back(eperp13(i,1));
    for (int i = 0; i < m.FN(); i++) vs.push_back(-eperp13(i,1));
    for (int i = 0; i < m.FN(); i++) vs.push_back(eperp21(i,1));
    for (int i = 0; i < m.FN(); i++) vs.push_back(-eperp21(i,1));
    for (int i = 0; i < m.FN(); i++) vs.push_back(eperp13(i,2));
    for (int i = 0; i < m.FN(); i++) vs.push_back(-eperp13(i,2));
    for (int i = 0; i < m.FN(); i++) vs.push_back(eperp21(i,2));
    for (int i = 0; i < m.FN(); i++) vs.push_back(-eperp21(i,2));

    // create sparse gradient operator matrix
    Eigen::SparseMatrix<double> G(3*m.FN(), m.VN());
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < (int) vs.size(); ++i) {
        triplets.push_back(Eigen::Triplet<double>(rs[i],cs[i],vs[i]));
    }
    G.setFromTriplets(triplets.begin(), triplets.end());

    //Eigen::SparseMatrix<double> Dx = G.block(0, 0, F.rows(), V.rows());
    //Eigen::SparseMatrix<double> Dy = G.block(F.rows(), 0, F.rows(), V.rows());
    //Eigen::SparseMatrix<double> Dz = G.block(2 * F.rows(), 0, F.rows(), V.rows());
    Eigen::SparseMatrix<double> Dx = G.block(0, 0, m.FN(), m.VN());
    Eigen::SparseMatrix<double> Dy = G.block(m.FN(), 0, m.FN(), m.VN());
    Eigen::SparseMatrix<double> Dz = G.block(2 * m.FN(), 0, m.FN(), m.VN());

    D1 = F1.col(0).asDiagonal() * Dx + F1.col(1).asDiagonal() * Dy + F1.col(2).asDiagonal() * Dz;
    D2 = F2.col(0).asDiagonal() * Dx + F2.col(1).asDiagonal() * Dy + F2.col(2).asDiagonal() * Dz;

    D1.makeCompressed();
    D2.makeCompressed();
}

void SLIM::UpdateJRW()
{
    Mesh& m = energy->m;

    MatrixXd Ji(m.FN(), 4);
    MatrixXd uv = X();

    Ji.col(0) = D1 * uv.col(0);
    Ji.col(1) = D2 * uv.col(0);
    Ji.col(2) = D1 * uv.col(1);
    Ji.col(3) = D2 * uv.col(1);

    for (auto& f : m.face) {
        // Compute jacobian matrix by combining the inverse of the map from T_e to 3D and the map from T_e to uv
        Point2d u10 = (f.V(1)->T().P() - f.V(0)->T().P());
        Point2d u20 = (f.V(2)->T().P() - f.V(0)->T().P());
        Eigen::Matrix2d fp;
        fp << u10[0], u20[0], u10[1], u20[1];

        Eigen::Matrix2d ji;
        int k = tri::Index(m, f);
        ji(0, 0) = Ji(k, 0);
        ji(0, 1) = Ji(k, 1);
        ji(1, 0) = Ji(k, 2);
        ji(1, 1) = Ji(k, 3);

        J[f] = fp * fm_inv[f];

        std::cout << "=============" << std::endl;
        std::cout << J[f] << std::endl;
        std::cout << "+++++++++++++" << std::endl;
        std::cout << ji << std::endl;
        std::cout << "=============" << std::endl;

        J[f] = ji;

        // Now find closest rotation for this face
        Eigen::Matrix2d U, V;
        Eigen::Vector2d sigma;
        Eigen::JacobiSVD<Eigen::Matrix2d> svd;
        svd.compute(J[f], Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd.matrixU(); V = svd.matrixV(); sigma = svd.singularValues();
        // Make sure U * V^{T} is a proper rotation (https://mathoverflow.net/questions/86539/closest-3d-rotation-matrix-in-the-frobenius-norm-sense)
        //bool reflected = false;
        //if (U.determinant() * V.determinant() < 0) {
        //    U.col(U.cols() - 1) *= -1;
        //    reflected = true;
        //}
        R[f] = U * V.transpose();
        if (R[f].determinant() < 0) {
            V.col(V.cols() - 1) *= -1;
            R[f] = U * V.transpose();
            V.col(V.cols() - 1) *= -1;
        }

        //if (reflected) U.col(U.cols() - 1) *= -1;

        // Update weights (eq 28 in the paper)
        Eigen::Vector2d Sw;
        for (int i = 0; i < 2; ++i) {
            double den = sigma[i] - 1;
            if (std::abs(den) < 1e8) {
                Sw[i] = 1;
            } else {
                Sw[i] = std::sqrt((sigma[i] - std::pow(sigma[i], 3)) / den);
            }
        }
        W[f] = U * Sw.asDiagonal() * U.transpose();
    }
}

Eigen::MatrixXd SLIM::ComputeDescentDirection()
{
    // Update jacobians, rotations and weights for each face
    UpdateJRW();

    // Find solution to the proxy energy
    Eigen::MatrixXd p_k = MinimizeProxyEnergy();

    // The descent direction is the one that points towards the proxy solution from the current point
    return p_k - X();
}

Eigen::MatrixXd SLIM::MinimizeProxyEnergy()
{
    // solves using least squares (so build the normal equations system)
    Eigen::SparseMatrix<double> A(4 * m.FN(), 2 * m.VN());
    BuildA(A);
    Eigen::SparseMatrix<double> At = A.transpose();
    At.makeCompressed();
    Eigen::SparseMatrix<double> id_m(At.rows(), At.rows());
    id_m.setIdentity();
    // add proximal penalty
    Eigen::SparseMatrix<double> L = At * diagAreaVector.asDiagonal() * A + lambda * id_m; //add also a proximal term
    L.makeCompressed();
    Eigen::VectorXd rhs;
    BuildRhs(At, rhs);
    L.makeCompressed();

    Eigen::VectorXd sol;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    sol = solver.compute(L).solve(rhs);

    assert((solver.info() == Eigen::Success) && "SLIM::MinimizeProxyEnergy solve failed");

    Eigen::MatrixXd solM(m.VN(), 2);
    solM.col(0) = sol.block(0, 0, m.VN(), 1);
    solM.col(1) = sol.block(m.VN(), 0, m.VN(), 1);
    return solM;
}

void SLIM::BuildA(Eigen::SparseMatrix<double>& A)
{
    // formula (35) in paper
    std::vector<Eigen::Triplet<double>> IJV;
    IJV.reserve(4 * (D1.outerSize() + D2.outerSize()));

    /*A = [W11*Dx, W12*Dx;
           W11*Dy, W12*Dy;
           W21*Dx, W22*Dx;
           W21*Dy, W22*Dy];*/
    for (int k = 0; k < D1.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(D1, k); it; ++it) {
            int dx_r = it.row();
            int dx_c = it.col();
            double val = it.value();

            const Eigen::Matrix2d& W_f = W[m.face[dx_r]];

            //IJV.push_back(Eigen::Triplet<double>(dx_r, dx_c, val * s.W_11(dx_r)));
            //IJV.push_back(Eigen::Triplet<double>(dx_r, s.v_n + dx_c, val * s.W_12(dx_r)));

            //IJV.push_back(Eigen::Triplet<double>(2 * s.f_n + dx_r, dx_c, val * s.W_21(dx_r)));
            //IJV.push_back(Eigen::Triplet<double>(2 * s.f_n + dx_r, s.v_n + dx_c, val * s.W_22(dx_r)));

            IJV.push_back(Eigen::Triplet<double>(dx_r, dx_c, val * W_f(0,0)));
            IJV.push_back(Eigen::Triplet<double>(dx_r, m.VN() + dx_c, val * W_f(0,1)));

            IJV.push_back(Eigen::Triplet<double>(2 * m.FN() + dx_r, dx_c, val * W_f(1,0)));
            IJV.push_back(Eigen::Triplet<double>(2 * m.FN() + dx_r, m.VN() + dx_c, val * W_f(1,1)));
        }
    }

    for (int k = 0; k < D2.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(D2, k); it; ++it) {
            int dy_r = it.row();
            int dy_c = it.col();
            double val = it.value();

            const Eigen::Matrix2d& W_f = W[m.face[dy_r]];

            //IJV.push_back(Eigen::Triplet<double>(s.f_n + dy_r, dy_c, val * s.W_11(dy_r)));
            //IJV.push_back(Eigen::Triplet<double>(s.f_n + dy_r, s.v_n + dy_c, val * s.W_12(dy_r)));

            //IJV.push_back(Eigen::Triplet<double>(3 * s.f_n + dy_r, dy_c, val * s.W_21(dy_r)));
            //IJV.push_back(Eigen::Triplet<double>(3 * s.f_n + dy_r, s.v_n + dy_c, val * s.W_22(dy_r)));

            IJV.push_back(Eigen::Triplet<double>(m.FN() + dy_r, dy_c, val * W_f(0,0)));
            IJV.push_back(Eigen::Triplet<double>(m.FN() + dy_r, m.VN() + dy_c, val * W_f(0,1)));

            IJV.push_back(Eigen::Triplet<double>(3 * m.VN() + dy_r, dy_c, val * W_f(1,0)));
            IJV.push_back(Eigen::Triplet<double>(3 * m.VN() + dy_r, m.VN() + dy_c, val * W_f(1,1)));
        }
    }

    A.setFromTriplets(IJV.begin(), IJV.end());
}

void SLIM::BuildRhs(const Eigen::SparseMatrix<double>& At, Eigen::VectorXd& rhs)
{
    Eigen::VectorXd f_rhs(4 * m.FN());
    f_rhs.setZero();
    /*b = [W11*R11 + W12*R21; (formula (36))
           W11*R12 + W12*R22;
           W21*R11 + W22*R21;
           W21*R12 + W22*R22];*/
    for (auto const& f : m.face) {
        int i = tri::Index(m, f);
        const Eigen::Matrix2d& W_f = W[f];
        const Eigen::Matrix2d& R_f = R[f];

        //f_rhs(i + 0 * m.FN()) = s.W_11(i) * s.Ri(i, 0) + s.W_12(i) * s.Ri(i, 1);
        //f_rhs(i + 1 * m.FN()) = s.W_11(i) * s.Ri(i, 2) + s.W_12(i) * s.Ri(i, 3);
        //f_rhs(i + 2 * m.FN()) = s.W_21(i) * s.Ri(i, 0) + s.W_22(i) * s.Ri(i, 1);
        //f_rhs(i + 3 * m.FN()) = s.W_21(i) * s.Ri(i, 2) + s.W_22(i) * s.Ri(i, 3);

        f_rhs(i + 0 * m.FN()) = W_f(0,0) * R_f(0,0) + W_f(0,1) * R_f(1,0);
        f_rhs(i + 1 * m.FN()) = W_f(0,0) * R_f(0,1) + W_f(0,1) * R_f(1,1);
        f_rhs(i + 2 * m.FN()) = W_f(1,0) * R_f(0,0) + W_f(1,1) * R_f(1,0);
        f_rhs(i + 3 * m.FN()) = W_f(1,0) * R_f(0,1) + W_f(1,1) * R_f(1,1);
    }
    Eigen::VectorXd uv_flat(2 * m.VN());
    for (auto const& v : m.vert) {
        int j = tri::Index(m, v);
        uv_flat(m.VN() * 0 + j) = v.T().U();
        uv_flat(m.VN() * 1 + j) = v.T().V();
    }

    rhs = (At * diagAreaVector.asDiagonal() * f_rhs + lambda * uv_flat);
}




