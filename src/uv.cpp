#include "uv.h"
#include "metric.h"
#include "math_utils.h"

using namespace vcg;

void StoreWedgeTexCoordAsAttribute(Mesh &m)
{
    //TODO assert MeshType has wedge tex coord
    Mesh::PerFaceAttributeHandle<TexCoordStorage> WTCSh
            = tri::Allocator<Mesh>::GetPerFaceAttribute<TexCoordStorage>(m, "WedgeTexCoordStorage");

    double scale;
    DistortionMetric::ComputeAreaScale(m, scale, ParameterizationGeometry::Model);
    scale = std::sqrt(1.0 / scale);
    int count = 0;

    for (auto &f : m.face) {
        if (DistortionMetric::AreaUV(f) == 0) {
            Point2d u10, u20;
            LocalIsometry(f.P(1) - f.P(0), f.P(2) - f.P(0), u10, u20);
            u10 *= scale; u20 *= scale;
            double u0_u = -1.0 - (rand() / (double) RAND_MAX);
            double u0_v = rand() / (double) RAND_MAX;
            Point2d u0{u0_u, u0_v};
            f.WT(0).P() = u0;
            f.WT(1).P() = u0 + u10;
            f.WT(2).P() = u0 + u20;
            /*WTCSh[&f].tc[0].P() = u0;
            WTCSh[&f].tc[0].N() = f.WT(0).N();
            WTCSh[&f].tc[1].P() = u0 + u10;
            WTCSh[&f].tc[1].N() = f.WT(1).N();
            WTCSh[&f].tc[2].P() = u0 + u20;
            WTCSh[&f].tc[2].N() = f.WT(2).N();*/
            assert((u10 ^ u20) > 0);
            count++;
        }// else {
            WTCSh[&f].tc[0].P() = f.WT(0).P();
            WTCSh[&f].tc[0].N() = f.WT(0).N();
            WTCSh[&f].tc[1].P() = f.WT(1).P();
            WTCSh[&f].tc[1].N() = f.WT(1).N();
            WTCSh[&f].tc[2].P() = f.WT(2).P();
            WTCSh[&f].tc[2].N() = f.WT(2).N();
       // }
    }
    // faces that were parameterized to degenerate triangles
    std::cout << "[LOG] Parameterized " << count << " zero uv area faces" << std::endl;
}

Box2d UVBox(const Mesh& m)
{
    Box2d uvbox;
    for(auto const& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            uvbox.Add(f.cWT(i).P());
        }
    }
    return uvbox;
}
