#include "arap.h"

#include "mesh_attribute.h"
#include "logging.h"
#include "math_utils.h"



ARAP::ARAP(Mesh& mesh)
    : m{mesh},
      max_iter{100}
{
}

void ARAP::FixVertex(Mesh::ConstVertexPointer vp, const vcg::Point2d& pos)
{
    fixed_i.push_back(tri::Index(m, vp));
    fixed_pos.push_back(pos);
}

void ARAP::FixBoundaryVertices()
{
    for (auto& v : m.vert) {
        if (v.IsB()) {
            fixed_i.push_back(tri::Index(m, v));
            fixed_pos.push_back(v.T().P());
        }
    }
}

int ARAP::FixSelectedVertices()
{
    int nfixed = 0;
    for (auto& v : m.vert) {
        if (v.IsS()) {
            fixed_i.push_back(tri::Index(m, v));
            fixed_pos.push_back(v.T().P());
            nfixed++;
        }
    }
    return nfixed;
}

/* This function fixes the vertices of an edge that is within 2pct of the target
 * edge length */
int ARAP::FixRandomEdgeWithinTolerance(double tol)
{
    std::unordered_set<int> fixed;
    for (int i : fixed_i)
        fixed.insert(i);

    auto tsa = GetTargetShapeAttribute(m);
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            double dcurr = (f.WT(i).P() - f.WT(f.Next(i)).P()).Norm();
            double dtarget = (tsa[f].P[i] - tsa[f].P[f.Next(i)]).Norm();
            if (std::abs((dcurr - dtarget) / dtarget) < tol) {
                if (fixed.count(tri::Index(m, f.V(i))) == 0 && fixed.count(tri::Index(m, f.V(f.Next(i)))) == 0) {
                    FixVertex(f.V(i), f.WT(i).P());
                    FixVertex(f.V(f.Next(i)), f.WT(f.Next(i)).P());
                    LOG_DEBUG << "Fixing vertices " << tri::Index(m, f.V(i)) << "   " << tri::Index(m, f.V(f.Next(i)));
                    return 2;
                }
            }
        }
    }
    return 0;
}

void ARAP::SetMaxIterations(int n)
{
    max_iter = n;
}

static std::vector<ARAP::Cot> ComputeCotangentVector(Mesh& m)
{
    std::vector<ARAP::Cot> cotan;
    cotan.reserve(m.FN());
    auto tsa = GetTargetShapeAttribute(m);
    double eps = std::numeric_limits<double>::epsilon();
    for (auto& f : m.face) {
        ARAP::Cot c;
        for (int i = 0; i < 3; ++i) {
            int j = (i+1)%3;
            int k = (i+2)%3;
            double alpha_i = std::max(VecAngle(tsa[f].P[j] - tsa[f].P[i], tsa[f].P[k] - tsa[f].P[i]), eps);
            c.v[i] = 0.5 * std::tan(M_PI_2 - alpha_i);
        }
        cotan.push_back(c);
    }
    return cotan;
}

void ARAP::ComputeSystemMatrix(Mesh& m, const std::vector<Cot>& cotan, Eigen::SparseMatrix<double>& L)
{
    using Td = Eigen::Triplet<double>;

    L.resize(m.VN(), m.VN());
    L.setZero();
    std::vector<Td> tri;
    auto Idx = [&m](const Mesh::VertexPointer vp) { return (int) tri::Index(m, vp); };
    for (auto &f : m.face) {
        int fi = tri::Index(m, f);
        for (int i = 0; i < 3; ++i) {
            if (std::find(fixed_i.begin(), fixed_i.end(), (int) tri::Index(m, f.V(i))) == fixed_i.end()) {
                Mesh::VertexPointer vi = f.V0(i);
                int j = (i+1)%3;
                Mesh::VertexPointer vj = f.V1(i);
                int k = (i+2)%3;
                Mesh::VertexPointer vk = f.V2(i);

                ensure_condition(Idx(vi) >= 0 && Idx(vi) < m.VN());
                ensure_condition(Idx(vj) >= 0 && Idx(vj) < m.VN());
                ensure_condition(Idx(vk) >= 0 && Idx(vk) < m.VN());

                double weight_ij = cotan[fi].v[k];
                double weight_ik = cotan[fi].v[j];

                if (!std::isfinite(weight_ij))
                    weight_ij = 1e-8;

                if (!std::isfinite(weight_ik))
                    weight_ik = 1e-8;

                tri.push_back(Td(Idx(vi), Idx(vj), -weight_ij));
                tri.push_back(Td(Idx(vi), Idx(vk), -weight_ik));
                tri.push_back(Td(Idx(vi), Idx(vi), (weight_ij + weight_ik)));
            }
        }
    }
    for (auto vi : fixed_i) {
        tri.push_back(Td(vi, vi, 1));
    }
    L.setFromTriplets(tri.begin(), tri.end());
    L.makeCompressed();
}

static std::vector<Eigen::Matrix2d> ComputeRotations(Mesh& m)
{
    auto tsa = GetTargetShapeAttribute(m);
    std::vector<Eigen::Matrix2d> rotations;
    rotations.reserve(m.FN());
    for (auto& f : m.face) {
        vcg::Point2d x10, x20;
        LocalIsometry(tsa[f].P[1] - tsa[f].P[0], tsa[f].P[2] - tsa[f].P[0], x10, x20);
        Eigen::Matrix2d Jf = ComputeTransformationMatrix(x10, x20, f.WT(1).P() - f.WT(0).P(), f.WT(2).P() - f.WT(0).P());
        Eigen::Matrix2d U, V;
        Eigen::Vector2d sigma;
        Eigen::JacobiSVD<Eigen::Matrix2d> svd;
        svd.compute(Jf, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd.matrixU(); V = svd.matrixV(); sigma = svd.singularValues();
        Eigen::MatrixXd R = U * V.transpose();
        if (R.determinant() < 0) {
            U.col(U.cols() - 1) *= -1;
            R = U * V.transpose();
        }
        rotations.push_back(R);
    }
    return rotations;
}

void ARAP::ComputeRHS(Mesh& m, const std::vector<Eigen::Matrix2d>& rotations, const std::vector<Cot>& cotan, Eigen::VectorXd& bu, Eigen::VectorXd& bv)
{
    auto Idx = [&m](const Mesh::VertexPointer vp) { return (int) tri::Index(m, vp); };
    bu = Eigen::VectorXd::Constant(m.VN(), 0);
    bv = Eigen::VectorXd::Constant(m.VN(), 0);
    auto tsa = GetTargetShapeAttribute(m);
    for (auto &f : m.face) {
        int fi = tri::Index(m, f);
        const Eigen::Matrix2d& Rf = rotations[fi];

        Eigen::Vector2d t[3];

        // TODO this should be computed once and stored in the object state
        Eigen::Vector2d x_10, x_20;
        LocalIsometry(tsa[f].P[1] - tsa[f].P[0], tsa[f].P[2] - tsa[f].P[0], x_10, x_20);
        t[0] = Eigen::Vector2d::Zero();
        t[1] = t[0] + x_10;
        t[2] = t[0] + x_20;

        for (int i = 0; i < 3; ++i) {
            Mesh::VertexPointer vi = f.V0(i);
            int j = (i+1)%3;
            int k = (i+2)%3;

            double weight_ij = cotan[fi].v[k];
            double weight_ik = cotan[fi].v[j];

            if (!std::isfinite(weight_ij))
                weight_ij = 1e-8;

            if (!std::isfinite(weight_ik))
                weight_ik = 1e-8;

            Eigen::Vector2d x_ij = t[i] - t[j];
            Eigen::Vector2d x_ik = t[i] - t[k];

            //Eigen::Vector2d x_ij(f.WT(i).U() - f.WT(j).U(), f.WT(i).V() - f.WT(j).V());
            //Eigen::Vector2d x_ik(f.WT(i).U() - f.WT(k).U(), f.WT(i).V() - f.WT(k).V());

            Eigen::Vector2d rhs = (weight_ij * Rf) * x_ij + (weight_ik * Rf) * x_ik;
            bu(Idx(vi)) += rhs.x();
            bv(Idx(vi)) += rhs.y();
        }
    }
    for (unsigned i = 0; i < fixed_i.size(); ++i) {
        bu(fixed_i[i]) = fixed_pos[i].X();
        bv(fixed_i[i]) = fixed_pos[i].Y();
    }
}

double ARAP::CurrentEnergy()
{
    double e = 0;
    double total_area = 0;
    auto tsa = GetTargetShapeAttribute(m);
    for (auto& f : m.face) {
        vcg::Point2d x10, x20;
        LocalIsometry(tsa[f].P[1] - tsa[f].P[0], tsa[f].P[2] - tsa[f].P[0], x10, x20);
        Eigen::Matrix2d Jf = ComputeTransformationMatrix(x10, x20, f.WT(1).P() - f.WT(0).P(), f.WT(2).P() - f.WT(0).P());
        Eigen::Matrix2d U, V;
        Eigen::Vector2d sigma;
        Eigen::JacobiSVD<Eigen::Matrix2d> svd;
        svd.compute(Jf, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd.matrixU(); V = svd.matrixV(); sigma = svd.singularValues();
        double area_f = 0.5 * ((tsa[f].P[1] - tsa[f].P[0]) ^ (tsa[f].P[2] - tsa[f].P[0])).Norm();
        total_area += area_f;
        e += area_f * (std::pow(sigma[0] - 1.0, 2.0) + std::pow(sigma[1] - 1.0, 2.0));
    }
    return e / total_area;
}

void ARAP::Solve()
{
    std::vector<Cot> cotan = ComputeCotangentVector(m);

    Eigen::SparseMatrix<double> A;
    ComputeSystemMatrix(m, cotan, A);

    Eigen::VectorXd xu = Eigen::VectorXd::Constant(m.VN(), 0);
    Eigen::VectorXd xv = Eigen::VectorXd::Constant(m.VN(), 0);

    double e = CurrentEnergy();

    // The system matrix is not symmetric
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    solver.analyzePattern(A);
    solver.factorize(A);

    if (solver.info() != Eigen::Success) {
        LOG_ERR << "Cotan matrix factorization failed: " << solver.info();
        assert(0);
    }

    /*
    std::string name = "arap_init.obj";
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            f.P(i) = Point3d(f.WT(i).U(), f.WT(i).V(), 0);
        }
    }
    tri::io::Exporter<Mesh>::Save(m, name.c_str(), tri::io::Mask::IOM_VERTCOORD);
    */

    bool converged = false;
    int iter = 0;
    while (!converged && iter < max_iter) {
        std::vector<Eigen::Matrix2d> rotations = ComputeRotations(m);
        Eigen::VectorXd bu(m.VN());
        Eigen::VectorXd bv(m.VN());
        ComputeRHS(m, rotations, cotan, bu, bv);

        Eigen::VectorXd xu_iter = solver.solve(bu);

        if (!(solver.info() == Eigen::Success)) {
            LOG_ERR << "ARAP solve failed";
            ensure_condition(0);
        }

        Eigen::VectorXd xv_iter = solver.solve(bv);

        if (!(solver.info() == Eigen::Success)) {
            LOG_ERR << "ARAP solve failed";
            ensure_condition(0);
        }

        for (auto& f : m.face) {
            for (int i = 0; i < 3; ++i) {
                int vi = tri::Index(m, f.V(i));
                f.WT(i).U() = xu_iter(vi);
                f.WT(i).V() = xv_iter(vi);
                f.V(i)->T() = f.WT(i);
            }
        }

        double e_curr = CurrentEnergy();

        double delta_e = e - e_curr;

        /*
        if (delta_e / e_curr < 1e-6) {
            LOG_DEBUG << "ARAP: convergence reached (change in the energy value is too small)";
            converged = true;
        }
        */

        xu = xu_iter;
        xv = xv_iter;
        e = e_curr;

        iter++;
    }

    if (iter == max_iter) {
        LOG_DEBUG << "ARAP: iteration limit reached";
    }

    // Extra step to ensure the fixed vertices do not move at all
    for (unsigned i = 0; i < fixed_i.size(); ++i) {
        m.vert[fixed_i[i]].T().P() = fixed_pos[i];
    }
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            f.WT(i) = f.cV(i)->T();
        }
    }

    LOG_DEBUG << "ARAP: Stopping after " << iter << " iterations, energy = " << e;
}


