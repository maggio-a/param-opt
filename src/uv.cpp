#include "uv.h"
#include "metric.h"
#include "math_utils.h"
#include "mesh_attribute.h"

using namespace vcg;
/*
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
            assert((u10 ^ u20) > 0);
            count++;
            */

void StoreWedgeTexCoordAsAttribute(Mesh &m)
{
    auto WTCSh = GetWedgeTexCoordStorageAttribute(m);
    for (auto &f : m.face) {
        WTCSh[&f].tc[0].P() = f.WT(0).P();
        WTCSh[&f].tc[0].N() = f.WT(0).N();
        WTCSh[&f].tc[1].P() = f.WT(1).P();
        WTCSh[&f].tc[1].N() = f.WT(1).N();
        WTCSh[&f].tc[2].P() = f.WT(2).P();
        WTCSh[&f].tc[2].N() = f.WT(2).N();
    }
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

double ComputeParameterizationScaleInfo(Mesh& m)
{
    auto info = GetParameterizationScaleInfoAttribute(m);
    for (auto& f : m.face) {
        double surfaceAreaFace = DistortionMetric::Area3D(f);
        double paramAreaFace = std::abs(DistortionMetric::AreaUV(f));
        if (surfaceAreaFace > 0.0 && paramAreaFace > 0.0) {
            info().surfaceArea += surfaceAreaFace;
            info().parameterArea += paramAreaFace;
            info().numNonZero++;
        }
    }
    info().scale = std::sqrt(info().parameterArea / info().surfaceArea);
}
