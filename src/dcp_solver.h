#ifndef DCPSOLVER_H
#define DCPSOLVER_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <type_traits>

#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>

#include "mesh.h"
#include "mesh_attribute.h"

#include <vcg/complex/complex.h>

/*
 * Discrete Conformal Parameterization
 * Reference: Desbrun, Meyer, Alliez - Intrinsic parameterizations of surface meshes
 * Parameterize a mesh by minimizing the Dirichlet's energy of the uv->3D map
 */
template <typename MeshType>
class DCPSolver
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

    bool Solve()
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
            /*
            tri::UpdateBounding<MeshType>::Box(mesh);
            const int bestAxis = mesh.bbox.MaxDim();
            for(VertexType &vv : mesh.vert) {
                if(vv.P()[bestAxis] <= mesh.bbox.min[bestAxis]) v0 = &vv;
                if(vv.P()[bestAxis] >= mesh.bbox.max[bestAxis]) v1 = &vv;
            }*/

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
                // as the matrix expression on the paper
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

#endif // DCPSOLVER_H

