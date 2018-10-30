#ifndef LINEAR_SOLVERS_H
#define LINEAR_SOLVERS_H

#include "mesh_attribute.h"

#include <vector>
#include <unordered_map>
#include <iostream>
#include <type_traits>

#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>

#include <vcg/complex/complex.h>


class LinearSolverInterface {
public:
    virtual ~LinearSolverInterface() { }
    virtual bool Solve() = 0;

};


/* Discrete Conformal Parameterizations
 * Reference: Desbrun, Meyer, Alliez - Intrinsic parameterizations of surface meshes
 * Parameterize a mesh by minimizing the Dirichlet's energy of the uv->3D map */
template <typename MeshType>
class DCPSolver : public LinearSolverInterface
{
public:

    using ScalarType = typename MeshType::ScalarType;
    using VertexType = typename MeshType::VertexType;
    using VertexPointer = typename MeshType::VertexPointer;
    using FaceType = typename MeshType::FaceType;
    using FacePointer = typename MeshType::FacePointer;
    using Coord3D = typename MeshType::CoordType;
    using CoordUV = typename FaceType::TexCoordType::PointType;

private:

    using Td = Eigen::Triplet<double>;

    MeshType &mesh;

    std::unordered_map<VertexPointer,CoordUV> constraints;
    std::vector<Td> coeffList;

    std::unordered_map<VertexPointer,int> vmap;

    int vn;

    int U(int i) { return i; }
    int V(int i) { return vn + i; }

public:

    DCPSolver(MeshType &m) : mesh{m}, vn{0} {}

    int U(VertexPointer v) { return U(vcg::tri::Index(mesh, v)); }
    int V(VertexPointer v) { return V(vcg::tri::Index(mesh, v)); }

    bool AddConstraint(VertexPointer vp, CoordUV uv)
    {
        return constraints.insert(std::make_pair(vp, uv)).second;
    }

    void BuildIndex()
    {
        for (auto& v : mesh.vert) {
            vmap[&v] = vn++;
        }
    }

    bool Solve() override
    {
        assert(HasTargetShapeAttribute(mesh));
        auto targetShape = GetTargetShapeAttribute(mesh);
        using namespace vcg;
        // Mark border vertices
        tri::UpdateTopology<MeshType>::FaceFace(mesh);
        tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(mesh);

        // Init index
        BuildIndex();

        if (constraints.size() < 2) {
            // Fix at least two vertices (needed since the energy is invariant to rotations
            // and translations, fixing two vertices acts as an anchor to the parametric space)
            VertexPointer v0 = nullptr;
            VertexPointer v1 = nullptr;

            vcg::Box3d box;
            for (auto& f : mesh.face) {
                for (int i = 0; i < 3; ++i) box.Add(VertPos(&f, i));
            }
            const int bestAxis = box.MaxDim();
            for (auto& f : mesh.face) {
                for (int i = 0; i < 3; ++i) {
                    Coord3D vpos = f.cP(i);
                    if (vpos[bestAxis] <= box.min[bestAxis]) v0 = f.V(i);
                    if (vpos[bestAxis] >= box.max[bestAxis]) v1 = f.V(i);
                }
            }
            assert((v0!=v1) && v0 && v1);
            AddConstraint(v0, CoordUV{0,0});
            AddConstraint(v1, CoordUV{1,1});
        }


        // Workaround to avoid inifinities if an angle is too small
        assert(std::is_floating_point<ScalarType>::value);
        ScalarType eps = std::numeric_limits<ScalarType>::epsilon();

        // Compute matrix coefficients
        coeffList.reserve(vn * 32);
        for (auto &f : mesh.face) {
            for (int i = 0; i < 3; ++i) if (constraints.find(f.V(i)) == constraints.end()) {
                // indices are given clockwise to avoid flipped parameterizations
                VertexPointer vi = f.V(i);
                Coord3D pi = targetShape[f][i];

                VertexPointer vj = f.V((i+2)%3);
                Coord3D pj = targetShape[f][(i+2)%3];

                VertexPointer vk = f.V((i+1)%3);
                Coord3D pk = targetShape[f][(i+1)%3];

                ScalarType alpha_ij = std::max(vcg::Angle(pi - pk, pj - pk), eps); // angle at k
                ScalarType alpha_ik = std::max(vcg::Angle(pi - pj, pk - pj), eps); // angle at j

                ScalarType weight_ij = std::tan(M_PI_2 - alpha_ij);
                ScalarType weight_ik = std::tan(M_PI_2 - alpha_ik);

                // Signs should be inverted wrt to the partial derivatives of the energy, but are the same
                // as the matrix expression in the paper
                // As long as the signs are consistent this should not change anything
                // since we solve for grad(f)(x) = 0, so its equivalent to multiply the system by -1
                coeffList.push_back(Td(U(vi), U(vj), weight_ij));
                coeffList.push_back(Td(U(vi), U(vk), weight_ik));
                coeffList.push_back(Td(U(vi), U(vi), -(weight_ij + weight_ik)));
                coeffList.push_back(Td(V(vi), V(vj), weight_ij));
                coeffList.push_back(Td(V(vi), V(vk), weight_ik));
                coeffList.push_back(Td(V(vi), V(vi), -(weight_ij + weight_ik)));

                if (vi->IsB()) {

                    // following from above, I invert the signs from eq (7) here too...
                    coeffList.push_back(Td(U(vi), V(vk), 1)); // pi/2 rotation: (x,y) -> (-y,x), hence V() as column index
                    coeffList.push_back(Td(U(vi), V(vj), -1));

                    coeffList.push_back(Td(V(vi), U(vk), -1));
                    coeffList.push_back(Td(V(vi), U(vj), 1));
                }
            }
        }

        // Prepare to solve

        // Allocate matrix and constants vector
        const int n = (vn + constraints.size()) * 2;
        Eigen::SparseMatrix<double,Eigen::ColMajor> A(n, n);
        Eigen::VectorXd b(Eigen::VectorXd::Zero(n));

        // Add fixed vertices coefficients
        int k = vn * 2;
        for (auto &c : constraints) {
            VertexPointer vp = c.first;
            CoordUV uv = c.second;
            coeffList.push_back(Td(U(vp), k, 1));
            coeffList.push_back(Td(k, U(vp), 1));
            b[k] = uv[0];
            k++;
            coeffList.push_back(Td(V(vp), k, 1));
            coeffList.push_back(Td(k, V(vp), 1));
            b[k] = uv[1];
            k++;
        }

        // Initialize matrix
        A.setFromTriplets(coeffList.begin(), coeffList.end());

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

        A.makeCompressed();
        solver.analyzePattern(A);
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            return false;
        }

        Eigen::VectorXd x = solver.solve(b);

        if (solver.info() != Eigen::Success) {
            return false;
        }

        // Copy back the texture coordinates
        vcg::Box2<typename CoordUV::ScalarType> uvBox;
        for (int i = 0; i < vn; ++i) {
            uvBox.Add(CoordUV(ScalarType(x[U(i)]), ScalarType(x[V(i)])));
        }

        double scale = 1.0 / std::max(uvBox.Dim().X(), uvBox.Dim().Y());

        for (auto &f : mesh.face) {
            for (int i = 0; i < 3; ++i) {
                VertexPointer vi = f.V(i);
                CoordUV uv(x[U(vi)], x[V(vi)]);
                f.WT(i).P() = (uv - uvBox.min) * scale;
                vi->T().P() = (uv - uvBox.min) * scale;
            }
        }

        return true;
    }

};


template <typename MeshType>
class UniformSolver : LinearSolverInterface
{
public:

    using ScalarType = typename MeshType::ScalarType;
    using VertexType = typename MeshType::VertexType;
    using VertexPointer = typename MeshType::VertexPointer;
    using FaceType = typename MeshType::FaceType;
    using FacePointer = typename MeshType::FacePointer;
    using Coord3D = typename MeshType::CoordType;
    using CoordUV = typename FaceType::TexCoordType::PointType;
    using IndexType = std::size_t;

private:

    using Td = Eigen::Triplet<double>;

    MeshType &mesh;

    std::unordered_map<IndexType,CoordUV> constraints;
    std::vector<Td> coeffList;

    int mode;

    int Idx(VertexPointer v)
    {
        return vcg::tri::Index(mesh, v);
    }

    bool AddConstraint(IndexType vi, CoordUV uv)
    {
        return constraints.insert(std::make_pair(vi, uv)).second;
    }

public:

    UniformSolver(MeshType &m) : mesh{m}, mode{0}
    {
        tri::UpdateTopology<MeshType>::FaceFace(m);
        assert(tri::Clean<MeshType>::MeshGenus(m) == 0);
    }

    void SetBoundaryMapProportional()
    {
        mode = 1;
    }

    void SetBoundaryMapUniform()
    {
        mode = 0;
    }

    bool Solve() override
    {
        // Mark border vertices
        tri::UpdateTopology<MeshType>::FaceFace(mesh);
        tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(mesh);


        std::vector<double> vTotalBorderLength;
        std::vector<std::vector<IndexType>> vBorderVertices;
        std::vector<std::vector<double>> vCumulativeBorder;

        tri::UpdateFlags<MeshType>::FaceClearV(mesh);
        for (auto& f : mesh.face) {
            for (int i = 0; i < 3; ++i) {
                if (!f.IsV() && face::IsBorder(f, i)) {
                    double totalBorderLength = 0;
                    std::vector<IndexType> borderVertices;
                    std::vector<double> cumulativeBorder;

                    face::Pos<FaceType> p(&f, i);
                    face::Pos<FaceType> startPos = p;
                    assert(p.IsBorder());
                    do {
                        assert(p.IsManifold());
                        p.F()->SetV();
                        borderVertices.push_back(tri::Index(mesh, p.V()));
                        cumulativeBorder.push_back(totalBorderLength);
                        totalBorderLength += EdgeLength(*p.F(), p.VInd());
                        p.NextB();
                    } while (p != startPos);
                    vTotalBorderLength.push_back(totalBorderLength);
                    vBorderVertices.push_back(borderVertices);
                    vCumulativeBorder.push_back(cumulativeBorder);
                }
            }
        }

        assert(vBorderVertices.size() == 1);
        // map border to the unit circle (store coord in vertex texcoord)
        constexpr double twoPi = 2 * M_PI;
        for (std::size_t i = 0; i < vBorderVertices[0].size(); ++i) {
            double angle;
            if (mode == 0) {
                angle = (i / double(vBorderVertices[0].size())) * twoPi;
            } else {
                assert(mode == 1) ;
                angle = (vCumulativeBorder[0][i] / vTotalBorderLength[0]) * twoPi;
            }

            Point2d uvCoord = Point2d{std::sin(angle), std::cos(angle)};
            mesh.vert[vBorderVertices[0][i]].T().P() = uvCoord;
            AddConstraint(vBorderVertices[0][i], uvCoord);
        }

        const int n = mesh.VN();
        Eigen::SparseMatrix<double,Eigen::ColMajor> A(n, n);
        Eigen::VectorXd b_u(Eigen::VectorXd::Zero(n));
        Eigen::VectorXd b_v(Eigen::VectorXd::Zero(n));

        std::unordered_map<VertexPointer, std::unordered_set<VertexPointer>> vadj;
        for (auto& f : mesh.face) {
            for (int i = 0; i < 3; ++i) {
                vadj[f.V0(i)].insert(f.V1(i));
                vadj[f.V0(i)].insert(f.V2(i));
            }
        }

        // Compute matrix coefficients
        coeffList.reserve(mesh.VN() * 32);
        for (auto& entry : vadj) {
            IndexType vi = tri::Index(mesh, entry.first);
            if (constraints.find(vi) == constraints.end()) {
                double c = entry.second.size();
                for (auto & vp : entry.second) {
                    coeffList.push_back(Td(Idx(entry.first), Idx(vp), - (1.0 / c)));
                }
            } else {
                b_u[vi] = constraints[vi][0];
                b_v[vi] = constraints[vi][1];
            }
            coeffList.push_back(Td(Idx(entry.first), Idx(entry.first), 1));
        }


        // Prepare to solve

        // Initialize matrix
        A.setFromTriplets(coeffList.begin(), coeffList.end());

        Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>, Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::StorageIndex>> solver;

        A.makeCompressed();
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            return false;
        }

        Eigen::VectorXd x_u = solver.solve(b_u);
        Eigen::VectorXd x_v = solver.solve(b_v);

        if (solver.info() != Eigen::Success) {
            return false;
        }

        // Copy back the texture coordinates
        vcg::Box2<typename CoordUV::ScalarType> uvBox;
        for (int i = 0; i < mesh.VN(); ++i) {
            uvBox.Add(CoordUV(x_u[i], x_v[i]));
        }

        double scale = 1.0 / std::max(uvBox.Dim().X(), uvBox.Dim().Y());

        for (auto &f : mesh.face) {
            for (int i = 0; i < 3; ++i) {
                VertexPointer vi = f.V(i);
                CoordUV uv(x_u[Idx(vi)], x_v[Idx(vi)]);
                f.WT(i).P() = (uv - uvBox.min) * scale;
                vi->T().P() = (uv - uvBox.min) * scale;
            }
        }

        return true;
    }

};


// TODO refactor using the mesh attributes...
template <typename MeshType>
class MeanValueSolver : LinearSolverInterface
{
public:

    using ScalarType = typename MeshType::ScalarType;
    using VertexType = typename MeshType::VertexType;
    using VertexPointer = typename MeshType::VertexPointer;
    using FaceType = typename MeshType::FaceType;
    using FacePointer = typename MeshType::FacePointer;
    using Coord3D = typename MeshType::CoordType;
    using CoordUV = typename FaceType::TexCoordType::PointType;
    using IndexType = std::size_t;

private:

    using Td = Eigen::Triplet<double>;

    MeshType &mesh;

    std::unordered_map<VertexPointer,CoordUV> constraints;
    std::vector<Td> coeffList;

    std::unordered_map<VertexPointer,int> vmap;

    int vn;

    int mode;

    int Idx(VertexPointer v)
    {
        return vcg::tri::Index(mesh, v);
    }

    bool AddConstraint(VertexPointer vp, CoordUV uv)
    {
        return constraints.insert(std::make_pair(vp, uv)).second;
    }

    void BuildIndex()
    {
        for (auto& v : mesh.vert) {
            vmap[&v] = vn++;
        }
    }

public:

    void UseCotangentWeights()
    {
        mode = 1;
    }

    void UseMeanValueWeights()
    {
        mode = 0;
    }

    MeanValueSolver(MeshType &m)
        : mesh{m},
          vn{0},
          mode{0}
    {
    }

    bool Solve() override
    {
        return Solve(DefaultVertexPosition<MeshType>{});
    }

private:

    template <typename VertexPosFct>
    bool Solve(VertexPosFct VertPos)
    {
        // Mark border vertices
        tri::UpdateTopology<MeshType>::FaceFace(mesh);
        tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(mesh);


        std::vector<double> vTotalBorderLength;
        std::vector<std::vector<IndexType>> vBorderVertices;
        std::vector<std::vector<double>> vCumulativeBorder;

        tri::UpdateFlags<MeshType>::FaceClearV(mesh);
        for (auto& f : mesh.face) {
            for (int i = 0; i < 3; ++i) {
                if (!f.IsV() && face::IsBorder(f, i)) {
                    double totalBorderLength = 0;
                    std::vector<IndexType> borderVertices;
                    std::vector<double> cumulativeBorder;

                    face::Pos<FaceType> p(&f, i);
                    face::Pos<FaceType> startPos = p;
                    assert(p.IsBorder());
                    do {
                        assert(p.IsManifold());
                        p.F()->SetV();
                        borderVertices.push_back(tri::Index(mesh, p.V()));
                        cumulativeBorder.push_back(totalBorderLength);
                        totalBorderLength += EdgeLength(*p.F(), p.VInd());
                        p.NextB();
                    } while (p != startPos);
                    vTotalBorderLength.push_back(totalBorderLength);
                    vBorderVertices.push_back(borderVertices);
                    vCumulativeBorder.push_back(cumulativeBorder);
                }
            }
        }

        assert(vBorderVertices.size() == 1);
        // map border to the unit circle (store coord in vertex texcoord)
        constexpr float twoPi = 2 * M_PI;
        for (std::size_t i = 0; i < vBorderVertices[0].size(); ++i) {
            //float angle = (vCumulativeBorder[0][i] / vTotalBorderLength[0]) * twoPi;
            float angle = (i / double(vBorderVertices[0].size())) * twoPi;
            Point2d uvCoord = Point2d{std::sin(angle), std::cos(angle)};
            mesh.vert[vBorderVertices[0][i]].T().P() = uvCoord;
            AddConstraint(&(mesh.vert[vBorderVertices[0][i]]), uvCoord);
        }

        const int n = mesh.VN();
        Eigen::SparseMatrix<double,Eigen::ColMajor> A(n, n);
        Eigen::VectorXd b_u(Eigen::VectorXd::Zero(n));
        Eigen::VectorXd b_v(Eigen::VectorXd::Zero(n));

        // Workaround to avoid inifinities if an angle is too small
        assert(std::is_floating_point<ScalarType>::value);
        ScalarType eps = std::numeric_limits<ScalarType>::epsilon();

        for (auto& f : mesh.face) {
            for (int i = 0; i < 3; ++i) if (constraints.find(f.V(i)) == constraints.end()) {
                VertexPointer vi = f.V(i);
                Coord3D pi = VertPos(&f, i);

                VertexPointer vj = f.V((i+2)%3);
                Coord3D pj = VertPos(&f, (i+2)%3);

                VertexPointer vk = f.V((i+1)%3);
                Coord3D pk = VertPos(&f, (i+1)%3);

                if (mode == 0) {
                    double angle_i = VecAngle(pi - pj, pi - pk);
                    if (!std::isfinite(angle_i)) {
                        angle_i = M_PI / 3.0;
                    }

                    double pij2 = (pi - pj).SquaredNorm();
                    double pik2 = (pi - pk).SquaredNorm();
                    double weight_ij = std::tan(angle_i / 2.0) / (pij2 > 0 ? pij2 : 1e-6)  ;
                    double weight_ik = std::tan(angle_i / 2.0) / (pik2 > 0 ? pik2 : 1e-6)  ;

                    if (!(std::isfinite(weight_ij) && std::isfinite(weight_ik))) {
                        return false;
                    }

                    coeffList.push_back(Td(Idx(vi), Idx(vj), weight_ij));
                    coeffList.push_back(Td(Idx(vi), Idx(vk), weight_ik));
                    coeffList.push_back(Td(Idx(vi), Idx(vi), -(weight_ij + weight_ik)));
                } else if (mode == 1) {
                    ScalarType alpha_ij = std::max(vcg::Angle(pi - pk, pj - pk), eps); // angle at k
                    ScalarType alpha_ik = std::max(vcg::Angle(pi - pj, pk - pj), eps); // angle at j

                    ScalarType weight_ij = std::tan(M_PI_2 - alpha_ij);
                    ScalarType weight_ik = std::tan(M_PI_2 - alpha_ik);

                    if (!(std::isfinite(weight_ij) && std::isfinite(weight_ik))) {
                        return false;
                    }

                    coeffList.push_back(Td(Idx(vi), Idx(vj), weight_ij));
                    coeffList.push_back(Td(Idx(vi), Idx(vk), weight_ik));
                    coeffList.push_back(Td(Idx(vi), Idx(vi), -(weight_ij + weight_ik)));

                } else {
                    assert(0);
                }
            }
        }

        // setup constraints coefficients
        for (auto& c : constraints) {
            VertexPointer vp = c.first;
            int vi = Idx(vp);
            coeffList.push_back(Td(vi, vi, 1));
            b_u[vi] = constraints[vp][0];
            b_v[vi] = constraints[vp][1];
        }

        A.setFromTriplets(coeffList.begin(), coeffList.end());
        Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>, Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::StorageIndex>> solver;

        A.makeCompressed();
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            return false;
        }

        Eigen::VectorXd x_u = solver.solve(b_u);
        Eigen::VectorXd x_v = solver.solve(b_v);

        if (solver.info() != Eigen::Success) {
            return false;
        }
        // Copy back the texture coordinates
        vcg::Box2<typename CoordUV::ScalarType> uvBox;
        for (int i = 0; i < mesh.VN(); ++i) {
            uvBox.Add(CoordUV(x_u[i], x_v[i]));
        }

        double scale = 1.0 / std::max(uvBox.Dim().X(), uvBox.Dim().Y());

        for (auto &f : mesh.face) {
            for (int i = 0; i < 3; ++i) {
                VertexPointer vi = f.V(i);
                CoordUV uv(x_u[Idx(vi)], x_v[Idx(vi)]);
                f.WT(i).P() = (uv - uvBox.min) * scale;
                vi->T().P() = (uv - uvBox.min) * scale;
            }
        }

        return true;
    }
};
#endif // LINEAR_SOLVERS_H

