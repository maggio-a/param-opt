/*
 * ATTRIBUTION
 * ===========
 *
 *
 * The code for Scalable locally injective mappings has been adapted from the
 * reference implementation that appears in libigl, available at
 *
 *   https://github.com/libigl/libigl
 *
 *
 *
 * The code for Composite majorization has been adapted from the implementation
 * publicly available at
 *
 *   https://github.com/Roipo/CompMajor
 *
 * */


#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "iterative_solvers.h"
#include "mesh.h"
#include "energy.h"
#include "math_utils.h"
#include "metric.h"
#include "timer.h"
#include "logging.h"


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
    //ensure_condition(x1 != 0 && x2 != 0);
    if (x2 < x1) std::swap(x1, x2);
    if (x1 > 0) return x1;
    else if (x2 > 0) return x2;
    else return std::numeric_limits<double>::max();
}


// Descent Method implementation
// =============================

DescentMethod::DescentMethod(std::shared_ptr<Energy> e)
    : m{e->m},
      energy{e}
{
    // empty
}

DescentMethod::~DescentMethod()
{
    // empty
}

bool DescentMethod::Iterate(double& gradientNorm, double& objValDiff, double& energyAfterIter)
{
    double energyPrev = energy->E();

    MatrixXd uv = X();
    MatrixXd grad = energy->Grad();

    MatrixXd dir;
    if (!ComputeDescentDirection(dir))
        return false;

    Search(uv, grad, dir);

    double energyCurr = energy->E();
    MatrixXd newGrad = energy->Grad();
    gradientNorm = std::sqrt(newGrad.cwiseProduct(newGrad).sum());
    objValDiff = energyPrev - energyCurr;

    energyAfterIter = energyCurr;

    return true;
}

void DescentMethod::UpdateCache()
{
    energy->UpdateCache();
}

double DescentMethod::SearchStrongWolfe(const MatrixXd& uv, const MatrixXd& grad, const MatrixXd& dir)
{
    double c1 = 0.2;
    double c2 = 0.9;
    double phi_zero = energy->E();
    double phiprime_zero = dir.cwiseProduct(grad).sum();

    double maxStepWithoutInversion = std::numeric_limits<double>::max();
    for (auto& f : m.face) {
        //ensure_condition(((f.V(1)->T().P() - f.V(0)->T().P()) ^ (f.V(2)->T().P() - f.V(0)->T().P())) > 0);
        double localStep = ComputeStepSizeNoFlip(uv.row(Index(m, f.V(0))), uv.row(Index(m, f.V(1))), uv.row(Index(m, f.V(2))),
                                                 dir.row(Index(m, f.V(0))), dir.row(Index(m, f.V(1))), dir.row(Index(m, f.V(2))));
        if (localStep < maxStepWithoutInversion) maxStepWithoutInversion = localStep;
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
            ensure_condition(alpha_curr < maxStepWithoutInversion);
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
                ensure_condition(alpha_curr < maxStepWithoutInversion);
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

    ensure_condition(alpha_curr < maxStepWithoutInversion);
    return alpha_curr;
}

double DescentMethod::Search(const MatrixXd& uv, const MatrixXd& grad, const MatrixXd& dir)
{
    double E_curr = energy->E();

    double c = 0.2;
    double descentCoeff = c * dir.cwiseProduct(grad).sum();

    //ensure_condition(descentCoeff < 0); // otherwise dir is not a descent direction
    if (descentCoeff >= 0)
        LOG_WARN << "Search direction is not a descent direction";

    double t = std::numeric_limits<double>::max();
    for (auto& f : m.face) {
        double tFace = ComputeStepSizeNoFlip(
                    uv.row(Index(m, f.V(0))), uv.row(Index(m, f.V(1))), uv.row(Index(m, f.V(2))),
                    dir.row(Index(m, f.V(0))), dir.row(Index(m, f.V(1))), dir.row(Index(m, f.V(2))));
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

    return E_t;
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
    ensure_condition(x.rows() == m.VN() && x.cols() == 2);
    for (auto& v : m.vert) {
        auto uv = x.row(Index(m, v));
        v.T().U() = uv(0);
        v.T().V() = uv(1);
    }
    /*
    for (auto& f : m.face) {
        Point2d u10 = f.V(1)->T().P() - f.V(0)->T().P();
        Point2d u20 = f.V(2)->T().P() - f.V(0)->T().P();
        ensure_condition((u10 ^ u20) > 0);
    }
    */
}


// Gradient Descent
// ================

GradientDescent::GradientDescent(std::shared_ptr<Energy> energy)
    : DescentMethod(energy)
{
    energy->CorrectScale();
}

bool GradientDescent::ComputeDescentDirection(Eigen::MatrixXd& dir)
{
    dir = - energy->Grad();
    return true;
}


// L-BFGS
// ======

LBFGS::LBFGS(std::shared_ptr<Energy> energy, std::size_t memory)
    : DescentMethod(energy),
      m{memory},
      sVec{},
      yVec{},
      rho{}
{
    energy->CorrectScale();
}

bool LBFGS::ComputeDescentDirection(Eigen::MatrixXd& dir)
{
    ensure_condition(sVec.size() == yVec.size());
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
    dir = -q;
    return true;
}

bool LBFGS::Iterate(double& gradientNorm, double& objValDiff, double& energyAfterIter)
{
    //energy->CorrectScale();
    double energyPrev = energy->E();
    MatrixXd uv = X();
    MatrixXd grad = energy->Grad();

    MatrixXd dir;
    if (!ComputeDescentDirection(dir))
        return false;

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

    energyAfterIter = energyCurr;
    return true;
}

/* Since the dimension of the problem changes, simply clear the history. Note
 * that this will have a rather large impact on the convergence rate of the
 * method if it is used too often */
void LBFGS::UpdateCache()
{
    DescentMethod::UpdateCache();
    sVec.clear();
    yVec.clear();
    rho.clear();
}


// Scalable locally injective mappings
// ===================================

SLIM::SLIM(std::shared_ptr<SymmetricDirichletEnergy> sd)
    : DescentMethod(sd),
      fm_inv{sd->m.face},
      J{sd->m.face},
      R{sd->m.face},
      W{sd->m.face},
      D1{},
      D2{},
      diagAreaVector{},
      lambda{0.0001},
      solver{},
      firstSolve{true},
      sd_energy{sd}
{
    UpdateCache();
}

bool SLIM::ComputeDescentDirection(Eigen::MatrixXd& dir)
{
    // Update jacobians, rotations and weights for each face
    UpdateJRW();

    // Find solution to the proxy energy minimization
    Eigen::MatrixXd p_k;
    if (!MinimizeProxyEnergy(p_k))
        return false;

    // The descent direction is the one that points towards the proxy solution from the current point
    for (auto const& v : m.vert) {
        p_k.row(tri::Index(m, v)) -= Eigen::Vector2d{v.T().U(), v.T().V()};
    }
    dir = p_k;
    return true;
}

void SLIM::UpdateCache()
{
    DescentMethod::UpdateCache();

    // Clear
    fm_inv.UpdateSize();
    J.UpdateSize();
    R.UpdateSize();
    W.UpdateSize();
    D1 = Eigen::SparseMatrix<double>{};
    D2 = Eigen::SparseMatrix<double>{};
    firstSolve = true;

    // Find 2d coords for the vertices of the faces relative to the frame that
    // is coplanar with the surface and has the origin in P0. x0 is (0,0), x1 is
    // (lenght(p1-p0), 0) and x2 is (1,0) rotated by the angle betwen (p1-p0)
    // and (p2-p0)
    for (auto& f : m.face) {
        Point3d p10 = energy->P1(&f) - energy->P0(&f);
        Point3d p20 = energy->P2(&f) - energy->P0(&f);
        Eigen::Matrix2d fm;
        Eigen::Vector2d x1, x2;
        LocalIsometry(p10, p20, x1, x2);
        fm.col(0) = x1;
        fm.col(1) = x2;
        // at this point the transformation from the canonical triangle (0,0) (1,0), (0,1) to the
        // mesh face is encoded in the matrix [x1 | x2], and we want the inverse
        fm_inv[f] = fm.inverse();
    }

    // Update rest of the solver state
    PrecomputeD12(m, D1, D2);
    diagAreaVector.resize(4 * m.FN());
    for (auto const& f : m.face) {
        double area;
        if (!f.IsScaffold())
            area = sd_energy->FaceArea(&f);
        else
            area = sd_energy->GetScaffoldFaceArea();
        int j = tri::Index(m, f);
        diagAreaVector(0 * m.FN() + j) = area;
        diagAreaVector(1 * m.FN() + j) = area;
        diagAreaVector(2 * m.FN() + j) = area;
        diagAreaVector(3 * m.FN() + j) = area;
    }
}

void SLIM::PrecomputeD12(const Mesh& m, Eigen::SparseMatrix<double>& D1, Eigen::SparseMatrix<double>& D2)
{
    Eigen::Matrix<double, Eigen::Dynamic, 3> eperp21(m.FN(), 3), eperp13(m.FN(), 3);
    Eigen::Matrix<double, Eigen::Dynamic, 3> F1(m.FN(), 3), F2(m.FN(), 3); // basis vectors for the tangent plane at each face

    //for (int i=0;i<F.rows();++i)
    for (auto const& f : m.face)
    {
        // #F x 3 matrices of triangle edge vectors, named after opposite vertices
        Eigen::Matrix<double, 1, 3> v32; (energy->P(&f, 2) - energy->P(&f, 1)).ToEigenVector(v32); // V.row(i3) - V.row(i2);
        Eigen::Matrix<double, 1, 3> v13; (energy->P(&f, 0) - energy->P(&f, 2)).ToEigenVector(v13); // V.row(i1) - V.row(i3);
        Eigen::Matrix<double, 1, 3> v21; (energy->P(&f, 1) - energy->P(&f, 0)).ToEigenVector(v21); // V.row(i2) - V.row(i1);
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
        eperp21.row(i) *= norm21 / dblA;
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
    triplets.reserve(vs.size());
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
    for (auto& f : m.face) {
        // Compute jacobian matrix by combining the inverse of the map from T_e to 3D and the map from T_e to uv
        Point2d u10 = (f.V(1)->T().P() - f.V(0)->T().P());
        Point2d u20 = (f.V(2)->T().P() - f.V(0)->T().P());
        Eigen::Matrix2d fp;
        fp << u10[0], u20[0], u10[1], u20[1];
        J[f] = fp * fm_inv[f];

        // Now find closest rotation for this face
        Eigen::Matrix2d U, V;
        Eigen::Vector2d sigma;
        Eigen::JacobiSVD<Eigen::Matrix2d> svd;
        svd.compute(J[f], Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd.matrixU(); V = svd.matrixV(); sigma = svd.singularValues();

        // Make sure U * V^T is a proper rotation (https://mathoverflow.net/questions/86539/closest-3d-rotation-matrix-in-the-frobenius-norm-sense)
        R[f] = U * V.transpose();
        if (R[f].determinant() < 0) {
            U.col(U.cols() - 1) *= -1;
            R[f] = U * V.transpose();
            U.col(U.cols() - 1) *= -1;
        }

        double weight = 1.0;
        if (f.IsScaffold())
            weight = sd_energy->GetScaffoldWeight();

        // Update weights (eq 28 in the paper)
        Eigen::Vector2d Sw;
        for (int i = 0; i < 2; ++i) {
            double den = sigma[i] - 1;
            if (std::abs(den) < 1e-8) {
                //Sw[i] = 1;
                Sw[i] = weight * 1;
            } else {
                //Sw[i] = std::sqrt((sigma[i] - std::pow(sigma[i], -3)) / den);
                Sw[i] = std::sqrt((weight * (sigma[i] - std::pow(sigma[i], -3))) / den);
            }
        }
        W[f] = U * Sw.asDiagonal() * U.transpose();
    }
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

            IJV.push_back(Eigen::Triplet<double>(m.FN() + dy_r, dy_c, val * W_f(0,0)));
            IJV.push_back(Eigen::Triplet<double>(m.FN() + dy_r, m.VN() + dy_c, val * W_f(0,1)));

            IJV.push_back(Eigen::Triplet<double>(3 * m.FN() + dy_r, dy_c, val * W_f(1,0)));
            IJV.push_back(Eigen::Triplet<double>(3 * m.FN() + dy_r, m.VN() + dy_c, val * W_f(1,1)));
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

#define U_IND(v) (tri::Index(m, v))
#define V_IND(v) (m.VN() + tri::Index(m, v))

#define MY_SET_U(v, val) \
        tripletVec.push_back(Td(U_IND(v), U_IND(v0), val*(-abac))); \
        tripletVec.push_back(Td(U_IND(v), U_IND(v1), val*(ab))); \
        tripletVec.push_back(Td(U_IND(v), U_IND(v2), val*(ac))); \
        tripletVec.push_back(Td(U_IND(v), V_IND(v0), val*(-ebec))); \
        tripletVec.push_back(Td(U_IND(v), V_IND(v1), val*(eb))); \
        tripletVec.push_back(Td(U_IND(v), V_IND(v2), val*(ec))); \
        rhs[U_IND(v)] += - (-ad-ef) * val;

#define MY_SET_V(v, val) \
        tripletVec.push_back(Td(V_IND(v), U_IND(v0), val*(-abac))); \
        tripletVec.push_back(Td(V_IND(v), U_IND(v1), val*(ab))); \
        tripletVec.push_back(Td(V_IND(v), U_IND(v2), val*(ac))); \
        tripletVec.push_back(Td(V_IND(v), V_IND(v0), val*(-ebec))); \
        tripletVec.push_back(Td(V_IND(v), V_IND(v1), val*(eb))); \
        tripletVec.push_back(Td(V_IND(v), V_IND(v2), val*(ec))); \
        rhs[V_IND(v)] += - (-ad-ef) * val;

inline void SetFromElementVertices(double a, double b, double c, double d, double e, double f, double area,
                                   const Mesh& m, Mesh::VertexPointer v0, Mesh::VertexPointer v1, Mesh::VertexPointer v2,
                                   std::vector<Eigen::Triplet<double>>& tripletVec, Eigen::VectorXd& rhs)
{
    using Td = Eigen::Triplet<double>;

    double ab = a*b;
    double ac = a*c;
    double ad = a*d;
    double eb = e*b;
    double ec = e*c;
    double ef = e*f;
    double abac = (ab + ac);
    double ebec = (eb + ec);

    MY_SET_U(v0, area*(-abac))
    MY_SET_U(v1, area*(ab))
    MY_SET_U(v2, area*(ac))
    MY_SET_V(v0, area*(-ebec))
    MY_SET_V(v1, area*(eb))
    MY_SET_V(v2, area*(ec))
}

void SLIM::BuildSystem(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &rhs)
{
    std::vector<Eigen::Triplet<double>> tripletVec(24 * m.FN());

    for (auto& f : m.face) {
        double area;
        if (!f.IsScaffold())
            area = sd_energy->FaceArea(&f);
        else
            area = sd_energy->GetScaffoldFaceArea();
        const Eigen::Matrix2d& Wf = W[f];
        const Eigen::Matrix2d& Qf = fm_inv[f];
        const Eigen::Matrix2d& Rf = R[f];
        SetFromElementVertices(Wf(0, 0), Qf(0, 0), Qf(1, 0), Rf(0, 0), Wf(0, 1), Rf(1, 0), area, m, f.V(0), f.V(1), f.V(2), tripletVec, rhs);
        SetFromElementVertices(Wf(0, 0), Qf(0, 1), Qf(1, 1), Rf(0, 1), Wf(0, 1), Rf(1, 1), area, m, f.V(0), f.V(1), f.V(2), tripletVec, rhs);
        SetFromElementVertices(Wf(1, 0), Qf(0, 0), Qf(1, 0), Rf(0, 0), Wf(1, 1), Rf(1, 0), area, m, f.V(0), f.V(1), f.V(2), tripletVec, rhs);
        SetFromElementVertices(Wf(1, 0), Qf(0, 1), Qf(1, 1), Rf(0, 1), Wf(1, 1), Rf(1, 1), area, m, f.V(0), f.V(1), f.V(2), tripletVec, rhs);
    }

    A.setZero();

    // TODO this is useless, I already know where the nonzero coefficients are in the matrix, so I can simply preallocate
    // space for those and set them directly
    A.setFromTriplets(tripletVec.begin(), tripletVec.end());
}

bool SLIM::MinimizeProxyEnergy(Eigen::MatrixXd& p_k)
{
    Timer t;
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

    Eigen::VectorXd sol;

    if (firstSolve) {
        solver.analyzePattern(L);
        firstSolve = false;
    }

    solver.factorize(L);
    if (!(solver.info() == Eigen::Success)) {
        LOG_ERR << "SLIM factorization failed";
        return false;
    }

    sol = solver.solve(rhs);
    if (!(solver.info() == Eigen::Success)) {
        LOG_ERR << "SLIM solve failed";
        return false;
    }

    p_k = Eigen::Map<MatrixXd>(sol.data(), m.VN(), 2);
    return true;
}


// Composite majorization
// ======================

CompMaj::CompMaj(std::shared_ptr<SymmetricDirichletEnergy> sd)
    : DescentMethod(sd),
      sd_energy(sd)
{
    UpdateCache();
}

void CompMaj::ComputeSurfaceGradientPerFace(Eigen::MatrixX3d& D1, Eigen::MatrixX3d& D2)
{
    //D1.resize(fn, 3);
    //D2.resize(fn, 3);

    Eigen::MatrixX3d F1(fn, 3);
    Eigen::MatrixX3d F2(fn, 3);

    Eigen::MatrixXd Dx(fn, 3), Dy(fn, 3), Dz(fn, 3);
    Eigen::Vector3i Pi;
    Pi << 1, 2, 0;
    Eigen::PermutationMatrix<3> P = Eigen::PermutationMatrix<3>(Pi);

    for (auto const& f : m.face)
    {
        int i = tri::Index(m, f);
        Eigen::Vector3d v32; (energy->P(&f, 2) - energy->P(&f, 1)).ToEigenVector(v32);
        Eigen::Vector3d v13; (energy->P(&f, 0) - energy->P(&f, 2)).ToEigenVector(v13);
        Eigen::Vector3d v21; (energy->P(&f, 1) - energy->P(&f, 0)).ToEigenVector(v21);
        Eigen::Vector3d v31; (energy->P(&f, 2) - energy->P(&f, 0)).ToEigenVector(v31);

        //Eigen::Vector3d n = v32.cross(-v13);
        Eigen::Vector3d n = v21.cross(v31);
        double dblA = std::sqrt(n.dot(n));
        Eigen::Vector3d u = n / dblA;

        F1.row(i) = v21.normalized();
        F2.row(i) = u.cross(v21).normalized();

        // #F x 3 matrices of triangle edge vectors, named after opposite vertices
        Eigen::Matrix3d e;
        e.col(0) = v21;
        e.col(1) = v32;
        e.col(2) = v13;

        Eigen::Matrix3d n_M;
        //n_M << 0, -n(2), n(1), n(2), 0, -n(0), -n(1), n(0), 0;
        n_M << 0, -u(2), u(1), u(2), 0, -u(0), -u(1), u(0), 0;
        Eigen::Matrix3d res = ((1. / dblA)*(n_M*e))*P;

        /* Dx, Dy and Dz are the components of the gradient vectors of each face at each vertex.
         * for face i, the gradient vector of vertex 0 is the triplet
         *   (Dx.row(i)(0), Dy.row(i)(0), Dz.row(i)(0))
         * and likewise for vertex 1 and vertex 2. For a geometric intuition of the
         * triangle gradients of the piecewise linear approximation of a scalar function
         * defined over a mesh cfr. Alec Jacobson's phd thesis
         *
         *   ``Algorithms and Interfaces for Real-Time Deformation of 2D and 3D Shapes''
         *
         * chapter 2.1 */
        Dx.row(i) = res.row(0);
        Dy.row(i) = res.row(1);
        Dz.row(i) = res.row(2);
    }
    /* The gradients are then transformed to the reference frame local to the triangle
     * (note that the gradient at each face is the linear combination of the gradient at
     * the three vertices */
    D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
    D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}

void CompMaj::UpdateCache()
{
    DescentMethod::UpdateCache();

    vn = m.VN();
    fn = m.FN();

    a.resize(fn);
    b.resize(fn);
    c.resize(fn);
    d.resize(fn);
    s.resize(fn, 2);
    v.resize(fn, 4);
    u.resize(fn, 4);

    Dsd[0].resize(6, fn);
    Dsd[1].resize(6, fn);

    Eigen::MatrixX3d D1, D2;
    ComputeSurfaceGradientPerFace(D1, D2);
    D1d = D1.transpose();
    D2d = D2.transpose();

    a1d.resize(6, fn);
    a2d.resize(6, fn);
    b1d.resize(6, fn);
    b2d.resize(6, fn);

    a1d.topRows(3) = 0.5 * D1d;
    a1d.bottomRows(3) = 0.5 * D2d;

    a2d.topRows(3) = -0.5 * D2d;
    a2d.bottomRows(3) = 0.5 * D1d;

    b1d.topRows(3) = 0.5 * D1d;
    b1d.bottomRows(3) = -0.5 * D2d;

    b2d.topRows(3) = 0.5 * D2d;
    b2d.bottomRows(3) = 0.5 * D1d;

    Hi.resize(fn);

    // prepare_hessian()

    II.clear();
    JJ.clear();
    auto PushPair = [&](int i, int j) {
        if (i > j)
            swap(i, j);
        II.push_back(i);
        JJ.push_back(j);
    };
    for (const auto& f : m.face) {
        // for every face there is a 6x6 local hessian
        // we only need the 21 values contained in the upper
        // triangle. they are access and also put into the
        // big hessian in column order.

        Eigen::Vector3i Fi;
        Fi << tri::Index(m, f.cV(0)), tri::Index(m, f.cV(1)), tri::Index(m, f.cV(2));
        // first column
        PushPair(Fi(0), Fi(0));

        // second column
        PushPair(Fi(0), Fi(1));
        PushPair(Fi(1), Fi(1));

        // third column
        PushPair(Fi(0), Fi(2));
        PushPair(Fi(1), Fi(2));
        PushPair(Fi(2), Fi(2));

        // fourth column
        PushPair(Fi(0), Fi(0) + vn);
        PushPair(Fi(1), Fi(0) + vn);
        PushPair(Fi(2), Fi(0) + vn);
        PushPair(Fi(0) + vn, Fi(0) + vn);

        // fifth column
        PushPair(Fi(0), Fi(1) + vn);
        PushPair(Fi(1), Fi(1) + vn);
        PushPair(Fi(2), Fi(1) + vn);
        PushPair(Fi(0) + vn, Fi(1) + vn);
        PushPair(Fi(1) + vn, Fi(1) + vn);

        // sixth column
        PushPair(Fi(0), Fi(2) + vn);
        PushPair(Fi(1), Fi(2) + vn);
        PushPair(Fi(2), Fi(2) + vn);
        PushPair(Fi(0) + vn, Fi(2) + vn);
        PushPair(Fi(1) + vn, Fi(2) + vn);
        PushPair(Fi(2) + vn, Fi(2) + vn);
    }
    SS = vector<double>(II.size(), 0.0);
}

void CompMaj::UpdateJ()
{
    for (auto& f : m.face) {
        int i = tri::Index(m, f);
        Eigen::Vector3d f_u, f_v;
        f_u << f.cV(0)->T().U(), f.cV(1)->T().U(), f.cV(2)->T().U();
        f_v << f.cV(0)->T().V(), f.cV(1)->T().V(), f.cV(2)->T().V();
        a(i) = D1d.col(i).transpose() * f_u;
        b(i) = D2d.col(i).transpose() * f_u;
        c(i) = D1d.col(i).transpose() * f_v;
        d(i) = D2d.col(i).transpose() * f_v;
    }
    detJuv = a.cwiseProduct(d) - b.cwiseProduct(c);
    for (int i = 0; i < fn; ++i) {
        MeshFace& f = m.face[i];
        if (detJuv[i] <= 0) {
            std::stringstream ss;
            ss << "Found negative determinant for ";
            if (f.IsMesh())
                ss << "mesh ";
            else if (f.IsHoleFilling())
                ss << " hole-filling ";
            else if (f.IsScaffold())
                ss << " scaffold ";
            else
                ensure_condition(0 && "impossible");
            ss << "face " << i << ", value = " << detJuv[i];
            LOG_ERR << ss.str();
        }
    }
    ensure_condition(!(detJuv.array() <= 0).any());
}

static inline void SSVD2x2(const Eigen::Matrix2d& A, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
    double e = (A(0) + A(3))*0.5;
    double f = (A(0) - A(3))*0.5;
    double g = (A(1) + A(2))*0.5;
    double h = (A(1) - A(2))*0.5;
    double q = std::sqrt((e*e) + (h*h));
    double r = std::sqrt((f*f) + (g*g));
    double a1 = std::atan2(g, f);
    double a2 = std::atan2(h, e);
    double rho = (a2 - a1)*0.5;
    double phi = (a2 + a1)*0.5;

    S(0) = q + r;
    S(1) = 0;
    S(2) = 0;
    S(3) = q - r;

    double c = std::cos(phi);
    double s = std::sin(phi);
    U(0) = c;
    U(1) = s;
    U(2) = -s;
    U(3) = c;

    c = std::cos(rho);
    s = std::sin(rho);
    V(0) = c;
    V(1) = -s;
    V(2) = s;
    V(3) = c;
}

void CompMaj::UpdateSSVDFunction()
{
    for (int i = 0; i < a.size(); ++i) {
        Eigen::Matrix2d J;
        J << a[i], b[i], c[i], d[i];
        Eigen::Matrix2d U, S, V;
        SSVD2x2(J, U, S, V);
        u.row(i) << U(0), U(1), U(2), U(3);
        v.row(i) << V(0), V(1), V(2), V(3);
        s.row(i) << S(0), S(3);
    }
}

void CompMaj::ComputeDenseSSVDDerivatives()
{
    //different columns belong to different faces
    Eigen::MatrixXd B(D1d * v.col(0).asDiagonal() + D2d * v.col(1).asDiagonal());
    Eigen::MatrixXd C(D1d * v.col(2).asDiagonal() + D2d * v.col(3).asDiagonal());

    Eigen::MatrixXd t1 = B * u.col(0).asDiagonal();
    Eigen::MatrixXd t2 = B * u.col(1).asDiagonal();
    Dsd[0].topRows(t1.rows()) = t1;
    Dsd[0].bottomRows(t1.rows()) = t2;
    t1 = C * u.col(2).asDiagonal();
    t2 = C * u.col(3).asDiagonal();
    Dsd[1].topRows(t1.rows()) = t1;
    Dsd[1].bottomRows(t1.rows()) = t2;
}

void CompMaj::ComputeHessian()
{
    // gradient of the sv energy formulation
    auto lambda1 = [](double val) { return 2.0 * val - 2.0 / (val*val*val); };
    // hessian of the sv energy formulation (it is a diagonal matrix)
    auto lambda2 = [](double val) { return 2.0 + 6.0 / (val*val*val*val); };

    Eigen::VectorXd areas(fn);
    for (int i = 0; i < areas.rows(); ++i) {
        MeshFace& f = m.face[i];
        if (!f.IsScaffold())
            areas(i) = sd_energy->FaceArea(&m.face[i]);
        else
            // I multiply here by the scaffold weight to simplify things
            areas(i) = sd_energy->GetScaffoldFaceArea() * sd_energy->GetScaffoldWeight();
    }

    Eigen::VectorXd gradS0 = areas.cwiseProduct(s.col(0).unaryExpr(lambda1));
    Eigen::VectorXd gradS1 = areas.cwiseProduct(s.col(1).unaryExpr(lambda1));
    Eigen::VectorXd HS0 = areas.cwiseProduct(s.col(0).unaryExpr(lambda2));
    Eigen::VectorXd HS1 = areas.cwiseProduct(s.col(1).unaryExpr(lambda2));

    Eigen::VectorXd alpha0 = 0.5 * (a + d);
    Eigen::VectorXd alpha1 = 0.5 * (c - b);
    Eigen::VectorXd beta0 = 0.5 * (a - d);
    Eigen::VectorXd beta1 = 0.5 * (c + b);

    for (const auto& f : m.face) {
        int i = tri::Index(m, f);
        Vector6d dS0i = Dsd[0].col(i);
        Vector6d dS1i = Dsd[1].col(i);
        Vector6d a1i = a1d.col(i);
        Vector6d a2i = a2d.col(i);
        Vector6d b1i = b1d.col(i);
        Vector6d b2i = b2d.col(i);
        /*
        Hi[i] = sd_energy->FaceArea(&f) * ComputeConvexConcaveFaceHessian(
                    a1i, a2i, b1i, b2i,
                    alpha0(i), alpha1(i), beta0(i), beta1(i),
                    dS0i, dS1i,
                    gradS0(i), gradS1(i),
                    HS0(i), HS1(i));
                    */
        Hi[i] = ComputeConvexConcaveFaceHessian(
                    a1i, a2i, b1i, b2i,
                    alpha0(i), alpha1(i), beta0(i), beta1(i),
                    dS0i, dS1i,
                    gradS0(i), gradS1(i),
                    HS0(i), HS1(i));

        // only use the values of the upper triangular portion of Hi
        int index = i * 21;
        for (int p = 0; p < 6; ++p) {
            for (int q = 0; q <= p; ++q) {
                SS[index++] = Hi[i](p, q);
            }
        }
    }

    for (int i = 0; i < fn; i++) {
        int base = 21 * i;
        SS[base] += 1e-6;
        SS[base + 2] += 1e-6;
        SS[base + 5] += 1e-6;
        SS[base + 9] += 1e-6;
        SS[base + 14] += 1e-6;
        SS[base + 20] += 1e-6;
    }
}

inline Matrix66d CompMaj::ComputeFaceConeHessian(const Vector6d& A1, const Vector6d& A2, double a1x, double a2x)
{
    double f2 = a1x * a1x + a2x * a2x;
    double invf =  1.0 / std::sqrt(f2);
    double invf3 = invf * invf * invf;

    Matrix66d A1A1t = A1 * A1.transpose();
    Matrix66d A2A2t = A2 * A2.transpose();
    Matrix66d A1A2t = A1 * A2.transpose();
    //Matrix66d A2A1t = A2 * A1.transpose();
    Matrix66d A2A1t = A1A2t.transpose();

    double a2 = a1x * a1x;
    double b2 = a2x * a2x;
    double ab = a1x * a2x;

    return (invf - invf3 * a2) * A1A1t + (invf - invf3 * b2) * A2A2t - invf3 * ab * (A1A2t + A2A1t);
}

inline Matrix66d CompMaj::ComputeConvexConcaveFaceHessian(
        const Vector6d& a1, const Vector6d& a2, const Vector6d& b1, const Vector6d& b2,
        double alpha0, double alpha1, double beta0, double beta1,
        const Vector6d& dS0i, const Vector6d& dS1i,
        double gradS0, double gradS1, double HS0, double HS1)
{
    Matrix66d H = HS0 * dS0i * dS0i.transpose() + HS1 * dS1i * dS1i.transpose();

    double walpha = gradS0 + gradS1;
    if (walpha > 0)
        H += walpha * ComputeFaceConeHessian(a1, a2, alpha0, alpha1);

    double wbeta = gradS0 - gradS1;
    if (wbeta > 1e-7)
        H += wbeta * ComputeFaceConeHessian(b1, b2, beta0, beta1);

    return H;
}

#include <wrap/io_trimesh/export.h>
bool CompMaj::ComputeDescentDirection(Eigen::MatrixXd& dir)
{
    UpdateJ();
    UpdateSSVDFunction();
    ComputeDenseSSVDDerivatives();
    ComputeHessian();

    {
        std::vector<Eigen::Triplet<double>> triplets;
        for (std::size_t i = 0; i < II.size(); ++i) {
            triplets.push_back(Eigen::Triplet<double>(II[i], JJ[i], SS[i]));
        }
        A.resize(2 * vn, 2 * vn);
        A.setFromTriplets(triplets.begin(), triplets.end());
    }

    Eigen::MatrixXd grad = energy->Grad();

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Upper> solver;

    solver.analyzePattern(A);

    solver.factorize(A);

    if (solver.info() != Eigen::Success) {
        LOG_ERR << "CompMaj factorization failed";
        tri::io::Exporter<Mesh>::Save(m, "factorization_failed.obj", tri::io::Mask::IOM_ALL);
        return false;
    }

    Eigen::VectorXd sol;
    sol = solver.solve(- Eigen::Map<Eigen::VectorXd>(grad.data(), 2 * vn));
    if (!(solver.info() == Eigen::Success)) {
        LOG_ERR << "CompMaj solve failed";
        return false;
    }

    dir = Eigen::Map<MatrixXd>(sol.data(), vn, 2);
    return true;
}

