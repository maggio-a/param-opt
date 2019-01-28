#include "cones.h"
#include "mesh.h"
#include "linear_solvers.h"
#include "metric.h"

#include <vector>


#include <wrap/io_trimesh/export.h>

int FindApproximateCone(Mesh& m, bool selectedFlag)
{

    std::vector<vcg::Point2d> savedVertexTexCoord(m.vert.size());
    std::vector<vcg::Point2d> savedWedgeTexCoord(m.face.size() * 3);

    for (unsigned i = 0; i < m.vert.size(); ++i)
        savedVertexTexCoord[i] = m.vert[i].T().P();
    for (unsigned i = 0; i < m.face.size(); ++i) {
        for (int k = 0; k < 3; ++k)
            savedWedgeTexCoord[i + k*m.face.size()] = m.face[i].WT(k).P();
    }


    // Compute a conformal parameterization

    DCPSolver<Mesh> solver(m);
    solver.Solve();


    // Detect the vertex index subject to the largest stretch

    double meshArea = 0;
    double parameterizationArea = 0;
    for (auto& f : m.face) {
        meshArea += DistortionMetric::Area3D(f);
        parameterizationArea += std::abs(DistortionMetric::AreaUV(f));
    }

    tri::UpdateTopology<Mesh>::VertexFace(m);
    std::vector<double> vertexScalings(m.vert.size(), 0.0);
    std::vector<int> vertexNumIncidentFaces(m.vert.size(), 0);




    std::vector<double> vert_3d_area(m.vert.size(), 0);
    std::vector<double> vert_uv_area(m.vert.size(), 0);



    for (auto& f : m.face) {
        double faceMeshAreaNormalized = DistortionMetric::Area3D(f) / meshArea;
        double faceParameterizationAreaNormalized = std::abs(DistortionMetric::AreaUV(f) / parameterizationArea);
        for (int i = 0; i < 3; ++i) {
            int vind = tri::Index(m, f.V(i));
            vertexScalings[vind] += (faceParameterizationAreaNormalized / faceMeshAreaNormalized);
            vertexNumIncidentFaces[vind]++;




            vert_3d_area[vind] += faceMeshAreaNormalized;
            vert_uv_area[vind] += faceParameterizationAreaNormalized;
        }
    }

    double minscale = std::numeric_limits<double>::max();
    int coneIndex = -1;
    for (unsigned i = 0; i < m.vert.size(); ++i) {
        /*
        double scale = (vertexScalings[i] / vertexNumIncidentFaces[i]);
        if (scale < minscale) {
            minscale = scale;
            coneIndex = i;
        }
        */

        if (!selectedFlag || m.vert[i].IsS()) {
            double scale2 = (vert_uv_area[i] / vert_3d_area[i]);

            if (scale2 < minscale) {
                minscale = scale2;
                coneIndex = i;
            }
        }
    }

    ensure_condition(coneIndex != -1);

    // Restore the original texture coordinates and return the vertex index

    for (unsigned i = 0; i < m.vert.size(); ++i)
        m.vert[i].T().P() = savedVertexTexCoord[i];
    for (unsigned i = 0; i < m.face.size(); ++i)
        for (int k = 0; k < 3; ++k)
            m.face[i].WT(k).P() = savedWedgeTexCoord[i + k*m.face.size()];

    return coneIndex;
}
