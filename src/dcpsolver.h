#ifndef DCPSOLVER_H
#define DCPSOLVER_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <type_traits>

#include <Eigen/Sparse>

#include "mesh.h"

/*
 * Discrete Conformal Parameterization
 * Reference: Desburn, Meyer, Alliez - Intrinsic parameterizations of surface meshes
 * Parameterize a mesh by minimizing the Dirichlet's energy of the uv->3D map
 */
template <typename MeshType>
class DCPSolver
{
    enum BorderManagement { skip, compute };

public:

    using ScalarType = typename MeshType::ScalarType;
    using VertexType = typename MeshType::VertexType;
    using VertexPointer = typename MeshType::VertexPointer;
    using FaceType = typename MeshType::FaceType;
    using FacePointer = typename MeshType::FacePointer;
    using Coord3D = typename MeshType::CoordType;
    using CoordUV = typename FaceType::TexCoordType::PointType;

    static const BorderManagement __trust_vertex_border_flags = skip;

private:

    using Td = Eigen::Triplet<double>;

    MeshType &mesh;

    std::unordered_map<VertexPointer,CoordUV> constraints;
    std::vector<Td> coeffList;

    std::unordered_map<VertexPointer,int> vmap;

    int vn;

    // support for DeclareUnique
    std::unordered_map<VertexPointer,VertexPointer> crossReferences;

    //int U(int i) { return 2*i; }
    //int V(int i) { return 2*i + 1; }
    int U(int i) { return i; }
    int V(int i) { return vn + i; }

public:

    DCPSolver(MeshType &m) : mesh{m}, vn{0} {}
   
    int U(VertexPointer v) { auto vi = crossReferences.find(v); return vi == crossReferences.end() ? U(vmap[v]) : U(vmap[vi->second]); }
    int V(VertexPointer v) { auto vi = crossReferences.find(v); return vi == crossReferences.end() ? V(vmap[v]) : V(vmap[vi->second]); }

    bool AddConstraint(VertexPointer vp, CoordUV uv)
    {
        return constraints.insert(std::make_pair(vp, uv)).second;
    }

    // function to map multiple vertices to the same index
    // no checks on duplicate entries, multiple declarations, etc...
    void DeclareUnique(const std::vector<VertexPointer>& vertices)
    {
        if (vertices.size() > 1) {
            for (auto vp : vertices) {
                crossReferences[vp] = vertices[0];
            }
        }
    }

    void BuildIndex()
    {
        for (auto& v : mesh.vert) {
            auto vi = crossReferences.find(&v);
            if (vi != crossReferences.end()) {
                VertexPointer rep = vi->second;
                if (vmap.count(rep) == 1) {
                    vmap[&v] = vmap[rep];
                }
                else {
                    vmap[&v] = vmap[rep] = vn++;
                }
            } else {
                vmap[&v] = vn++;
            }
        }
    }

    // Mock PoissonSolver interface
    void Init() {}
    void FixDefaultVertices() {}
    bool IsFeasible() { return true; }
    bool SolvePoisson() { return Solve(); }

    bool Solve(BorderManagement mode = compute)
    {
        // Mark border vertices
        if (mode != __trust_vertex_border_flags) {
            tri::UpdateTopology<MeshType>::FaceFace(mesh);
            tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(mesh);
        }
        // Init index
        BuildIndex();

        /*
        std::cout << "Border indices" << std::endl;
        for (auto p : vmap) {
            if (p.first->IsB()) std::cout << p.second << std::endl;
        }
        */

        // Workaround to avoid inifinities if an angle is too small
        assert(std::is_floating_point<ScalarType>::value);
        ScalarType eps = std::numeric_limits<ScalarType>::epsilon();

        // Compute matrix coefficients
        coeffList.reserve(vn * 32);
        for (auto &f : mesh.face) {
            for (int i = 0; i < 3; ++i) {
                VertexPointer vi = f.V0(i);
                VertexPointer vj = f.V1(i);
                VertexPointer vk = f.V2(i);

                ScalarType alpha_ij = std::max(vcg::Angle(vi->P() - vk->P(), vj->P() - vk->P()), eps);
                ScalarType weight_ij = ScalarType(1) / std::tan(alpha_ij);

                ScalarType alpha_ik = std::max(vcg::Angle(vi->P() - vj->P(), vk->P() - vj->P()), eps);
                ScalarType weight_ik = ScalarType(1) / std::tan(alpha_ik);

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
                    coeffList.push_back(Td(U(vi), V(vk), -1)); // pi/2 rotation: (x,y) -> (-y,x), hence V() as column index
                    coeffList.push_back(Td(U(vi), V(vj), 1));

                    coeffList.push_back(Td(V(vi), U(vk), -1));
                    coeffList.push_back(Td(V(vi), U(vj), 1));
                }
            }
        }

        // Fix at least two vertices (needed since the energy is invariant to rotations
        // and translations, fixing two vertices acts as an anchor to the parametric space)
        
        tri::UpdateBounding<MeshType>::Box(mesh);

        VertexPointer v0 = nullptr;
        VertexPointer v1 = nullptr;
        const int bestAxis = mesh.bbox.MaxDim();
        for(VertexType &vv : mesh.vert) {
            if(vv.P()[bestAxis] <= mesh.bbox.min[bestAxis]) v0 = &vv;
            if(vv.P()[bestAxis] >= mesh.bbox.max[bestAxis]) v1 = &vv;
        }
        assert( (v0!=v1) && v0 && v1);

        AddConstraint(v0, CoordUV{0,0});
        AddConstraint(v1, CoordUV{1,1});

        // Prepare to solve
        int n = (vn + constraints.size()) * 2;

        // Allocate matrix and constants vector
        Eigen::SparseMatrix<double> A(n, n);
        Eigen::VectorXd b(Eigen::VectorXd::Zero(n));

        // Add constraints to the system
        int h = vn * 2; // first constraint index
        for (auto &c : constraints) {
            VertexPointer vp = c.first;
            CoordUV uv = c.second;

            // constrain u
            coeffList.push_back(Td(U(vp), h, 1));
            coeffList.push_back(Td(h, U(vp), 1));
            b[h] = uv[0];
            h++;

            // constrain v
            coeffList.push_back(Td(V(vp), h, 1));
            coeffList.push_back(Td(h, V(vp), 1));
            b[h] = uv[1];
            h++;
        }

        // Initialize matrix
        A.setFromTriplets(coeffList.begin(), coeffList.end());

        //std::cout << Eigen::MatrixXd(A) << std::endl;
        //std::cout << b << std::endl;

        // Solve TODO find out why the conjugate gradient fails with x = NaNs
        //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholesky;

        cholesky.compute(A);
        if (cholesky.info() != Eigen::Success) {
            std::cout << "Factorization failed" << std::endl;
            return false;
        }

        Eigen::VectorXd x = cholesky.solve(b);

        //std::cout << "# iterations    " << cg.iterations() << std::endl;
        //std::cout << "estimated error " << cg.error() << std::endl;

        if (cholesky.info() != Eigen::Success) {
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

