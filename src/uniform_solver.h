#ifndef UNIFORM_SOLVER_H
#define UNIFORM_SOLVER_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <type_traits>

#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>

#include "mesh.h"
#include "metric.h"
#include "math_utils.h"

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clean.h>

template <typename MeshType>
class UniformSolver
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

public:

    UniformSolver(MeshType &m) : mesh{m}
    {
        tri::UpdateTopology<MeshType>::FaceFace(m);
        assert(tri::Clean<MeshType>::MeshGenus(m) == 0);
    }

private:

    int Idx(VertexPointer v) { return vcg::tri::Index(mesh, v); }

    bool AddConstraint(IndexType vi, CoordUV uv)
    {
        return constraints.insert(std::make_pair(vi, uv)).second;
    }

public:

    bool Solve()
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
            //float angle = (vCumulativeBorder[0][i] / vTotalBorderLength[0]) * twoPi;
            double angle = (i / double(vBorderVertices[0].size())) * twoPi;
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
                    //coeffList.push_back(Td(V(entry.first), V(vp), coeff));
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

        //std::cout << Eigen::MatrixXd(A) << std::endl;
        //std::cout << b << std::endl;

        //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver; requires to reorder vertices so boundary vs are last
        Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>, Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::StorageIndex>> solver;

        A.makeCompressed();
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            std::cout << "Factorization failed" << std::endl;
            return false;
        }

        Eigen::VectorXd x_u = solver.solve(b_u);
        Eigen::VectorXd x_v = solver.solve(b_v);

        if (solver.info() != Eigen::Success) {
            std::cout << "Solving failed" << std::endl;
            return false;
        }

        // Copy back the texture coordinates
        vcg::Box2<typename CoordUV::ScalarType> uvBox;
        for (int i = 0; i < mesh.VN(); ++i) {
            uvBox.Add(CoordUV(x_u[i], x_v[i]));
        }

        double scale = 1.0f / std::max(uvBox.Dim().X(), uvBox.Dim().Y());

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
                //assert(DistortionMetric::AreaUV(f) > 0);
            }
        }

        return true;
    }

};

#endif // UNIFORM_SOLVER_H

