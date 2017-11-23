#ifndef FIXEX_BORDER_BIJECTIVE_H
#define FIXEX_BORDER_BIJECTIVE_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <type_traits>

#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>

#include "mesh.h"

#include <vcg/complex/complex.h>

template <typename FaceType>
float EdgeLength(const FaceType& f, int i) {
    return vcg::Distance(f.cV0(i)->P(), f.cV1(i)->P());
}

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

    UniformSolver(MeshType &m) : mesh{m}, vn{0} {}

private:

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

public:

    bool Solve()
    {
        // Mark border vertices
        tri::UpdateTopology<MeshType>::FaceFace(mesh);
        tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(mesh);

        // Init index
        BuildIndex();

        float totalBorderLength = 0.0f;

        std::vector<VertexPointer> borderVertices;
        std::vector<float> cumulativeBorder;

        tri::UpdateFlags<MeshType>::FaceClearV(mesh);
        bool borderFound = false;
        for (auto& f : mesh.face) {
            for (int i = 0; i < 3; ++i) {
                if (!f.IsV() && face::IsBorder(f, i)) {
                    assert(!borderFound && "Multiple borders detected");
                    borderFound = true;
                    face::Pos<FaceType> p(&f, i);
                    face::Pos<FaceType> startPos = p;
                    assert(p.IsBorder());
                    do {
                        assert(p.IsManifold());
                        p.F()->SetV();
                        borderVertices.push_back(p.V());
                        cumulativeBorder.push_back(totalBorderLength);
                        totalBorderLength += EdgeLength(*p.F(), p.VInd());
                        p.NextB();
                    } while (p != startPos);
                }
            }
        }

        // map border to the unit circle (store coord in vertex texcoord)
        constexpr float twoPi = 2 * M_PI;
        for (std::size_t i = 0; i < borderVertices.size(); ++i) {
            float angle = (cumulativeBorder[i] / totalBorderLength) * twoPi;
            borderVertices[i]->T().P() = Point2f{std::sin(angle), std::cos(angle)};
            AddConstraint(borderVertices[i], borderVertices[i]->T().P());
        }

        // Workaround to avoid inifinities if an angle is too small
        assert(std::is_floating_point<ScalarType>::value);
        ScalarType eps = std::numeric_limits<ScalarType>::epsilon();

        const int n = vn * 2;
        Eigen::SparseMatrix<double,Eigen::ColMajor> A(n, n);
        Eigen::VectorXd b(Eigen::VectorXd::Zero(n));

        std::unordered_map<VertexPointer, std::unordered_set<VertexPointer>> vadj;
        for (auto& f : mesh.face) {
            for (int i = 0; i < 3; ++i) {
                vadj[f.V0(i)].insert(f.V1(i));
                vadj[f.V0(i)].insert(f.V2(i));
            }
        }

        // Compute matrix coefficients
        coeffList.reserve(vn * 32);
        for (auto& entry : vadj) {
            if (constraints.find(entry.first) == constraints.end()) {
                for (auto & vp : entry.second) {
                    float coeff = - 1.0f / entry.second.size();
                    coeffList.push_back(Td(U(entry.first), U(vp), coeff));
                    coeffList.push_back(Td(V(entry.first), V(vp), coeff));
                }
            } else {
                b[U(entry.first)] = constraints[entry.first][0];
                b[V(entry.first)] = constraints[entry.first][1];
            }
            coeffList.push_back(Td(U(entry.first), U(entry.first), 1));
            coeffList.push_back(Td(V(entry.first), V(entry.first), 1));
        }

        // Prepare to solve

        // Initialize matrix
        A.setFromTriplets(coeffList.begin(), coeffList.end());

        //std::cout << Eigen::MatrixXd(A) << std::endl;
        //std::cout << b << std::endl;

        Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>, Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::StorageIndex>> solver;

        A.makeCompressed();
        solver.analyzePattern(A);

        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            std::cout << "Factorization failed" << std::endl;
            return false;
        }

        Eigen::VectorXd x = solver.solve(b);

        if (solver.info() != Eigen::Success) {
            std::cout << "Solving failed" << std::endl;
            return false;
        }

        // Copy back the texture coordinates
        vcg::Box2<typename CoordUV::ScalarType> uvBox;
        for (int i = 0; i < vn; ++i) {
            uvBox.Add(CoordUV(ScalarType(x[U(i)]), ScalarType(x[V(i)])));
        }

        float scale = 1.0f / std::max(uvBox.Dim().X(), uvBox.Dim().Y());

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

#endif // DCPSOLVER_H

