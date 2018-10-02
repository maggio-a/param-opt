#ifndef MEAN_VALUE_PARAM_H
#define MEAN_VALUE_PARAM_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <type_traits>

#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>

#include "mesh.h"
#include "math_utils.h"
#include "vertex_position.h"

#include <vcg/complex/complex.h>

template <typename MeshType>
class MeanValueSolver
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

    int mode = 0;

public:

    void UseCotangentWeights()
    {
        mode = 1;
    }

    MeanValueSolver(MeshType &m) : mesh{m}, vn{0} {}

    int Idx(VertexPointer v) { return vcg::tri::Index(mesh, v); }

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

    bool Solve(float lambda = 1.0f, float mi = 0.0f)
    {
        return Solve(DefaultVertexPosition<MeshType>{});
    }

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
                        std::cout << "Failed to compute matrix coefficients for face " << tri::Index(mesh, f) << std::endl;
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
                        std::cout << "Failed to compute matrix coefficients for face " << tri::Index(mesh, f) << std::endl;
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
            std::cout << "Matrix factorization failed" << std::endl;
            return false;
        }

        Eigen::VectorXd x_u = solver.solve(b_u);
        Eigen::VectorXd x_v = solver.solve(b_v);

        // Copy back the texture coordinates
        vcg::Box2<typename CoordUV::ScalarType> uvBox;
        for (int i = 0; i < mesh.VN(); ++i) {
            uvBox.Add(CoordUV(x_u[i], x_v[i]));
        }

        double scale = 1.0 / std::max(uvBox.Dim().X(), uvBox.Dim().Y());

        for (auto &f : mesh.face) {
            for (int i = 0; i < 3; ++i) {
                VertexPointer vi = f.V(i);
                //CoordUV uv(x[Idx(vi)], x[V(vi)]);
                CoordUV uv(x_u[Idx(vi)], x_v[Idx(vi)]);
                f.WT(i).P() = (uv - uvBox.min) * scale;
                vi->T().P() = (uv - uvBox.min) * scale;
            }
        }

        for (auto &f : mesh.face) {
            if (DistortionMetric::AreaUV(f) <= 0) {
                std::cout << "Failed to compute injective parameterization" << std::endl;
                return false;
            }
        }

        return true;
    }
};

#endif // MEAN_VALUE_PARAM_H

