#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "iterative.h"
#include "mesh.h"
#include "energy.h"

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


double get_smallest_pos_quad_zero(double a,double b, double c) {
  double t1,t2;
  if (a != 0) {
    double delta_in = pow(b,2) - 4*a*c;
    if (delta_in < 0) {
      return INFINITY;
    }
    double delta = sqrt(delta_in);
    t1 = (-b + delta)/ (2*a);
    t2 = (-b - delta)/ (2*a);
  } else {
    t1 = t2 = -b/c;
  }
  assert (std::isfinite(t1));
  assert (std::isfinite(t2));

  double tmp_n = min(t1,t2);
  t1 = max(t1,t2); t2 = tmp_n;
  if (t1 == t2) {
    return INFINITY; // means the orientation flips twice = doesn't flip?
  }
  // return the smallest negative root if it exists, otherwise return infinity
  if (t1 > 0) {
    if (t2 > 0) {
      return t2;
    } else {
      return t1;
    }
  } else {
    return INFINITY;
  }
}

// this functions takes the 2d base coordinates and the offsets for each vertex, and returns the smallest positive
// step that prevents a triangle flip, or numeric_limits::max() if the step is unbounded
// Negli appunti pij = (u_ij, v_ij), uguale per pki
static double ComputeStepSizeNoFlip(const Vector2d& pi, const Vector2d& pj, const Vector2d& pk, const Vector2d& di, const Vector2d& dj, const Vector2d& dk)
{
//    Vector2d dji = dj - di; Vector2d dki = dk - di;
//    Vector2d pji = pj - pi; Vector2d pki = pk - pi;

//    double a = dji[0]*dki[1] - dji[1]*dki[0];
//    double b = pji[0]*dki[1] + pki[1]*dji[0] - pji[1]*dki[0] - pki[0]*dji[1];
//    double c = pji[0]*pki[1] - pji[1]*pki[0];

#define U11 pi[0]
#define U12 pi[1]
#define U21 pj[0]
#define U22 pj[1]
#define U31 pk[0]
#define U32 pk[1]

#define V11 di[0]
#define V12 di[1]
#define V21 dj[0]
#define V22 dj[1]
#define V31 dk[0]
#define V32 dk[1]


double a = V11*V22 - V12*V21 - V11*V32 + V12*V31 + V21*V32 - V22*V31;
double b = U11*V22 - U12*V21 - U21*V12 + U22*V11 - U11*V32 + U12*V31 + U31*V12 - U32*V11 + U21*V32 - U22*V31 - U31*V22 + U32*V21;
double c = U11*U22 - U12*U21 - U11*U32 + U12*U31 + U21*U32 - U22*U31;

return get_smallest_pos_quad_zero(a, b, c);


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

    double gamma = 0.8;
    double E_t = E_curr;
    int numIter = 0;

    t /= gamma;
    while (E_t > E_curr + t * descentCoeff) {
        t *= gamma;
        SetX(uv + (t * dir));
        E_t = energy->E();
        numIter++;
    }

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




