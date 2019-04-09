#include "mesh.h"
#include "mesh_graph.h"
#include "uv.h"
#include "texture_optimization.h"
#include "texture_rendering.h"
#include "timer.h"
#include "gl_utils.h"
#include "texture.h"

#include "logging.h"

#include <wrap/io_trimesh/io_mask.h>
#include <wrap/io_trimesh/export.h>

#include <vcg/complex/algorithms/update/color.h>

#include <string>
#include <vector>
#include <iostream>

#include <QCoreApplication>
#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>


using namespace vcg;


int main(int argc, char *argv[])
{
    // Make sure the executable directory is added to Qt's library path
    {
        QFileInfo executableInfo(argv[0]);
        QCoreApplication::addLibraryPath(executableInfo.dir().absolutePath());
    }
    
    GLInit();

    Mesh m;
    TextureObjectHandle textureObject;
    int loadMask;

    ensure_condition(argc > 1);

    if (LoadMesh(m, argv[1], textureObject, loadMask) == false) {
        std::cout << "Failed to open mesh " << argv[1] << std::endl;
        std::exit(-1);
    }

    ensure_condition(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    tri::UpdateTopology<Mesh>::FaceFace(m);

    std::string name = m.name;
    int vn = m.VN();
    int fn = m.FN();

    int nonManifoldEdgeCount = tri::Clean<Mesh>::CountNonManifoldEdgeFF(m);
    int nonManifoldVertexCount = tri::Clean<Mesh>::CountNonManifoldVertexFF(m);

    int textureNum = textureObject->ArraySize();

    auto graph = ComputeParameterizationGraph(m, textureObject);

    std::size_t atlasRegionCount = graph->Count();

    std::size_t emptyRegionCount = 0;
    for (auto& entry : graph->charts)
        if (entry.second->AreaUV() == 0.0)
            emptyRegionCount++;

    std::vector<RasterizedParameterizationStats> textureStats = GetRasterizationStats(m, textureObject);
    while (textureStats.size() < 16) {
        textureStats.push_back({});
    }

    std::cout << "NAME , VN , FN , NMANIF_EDGE , NMANIF_VERT , CHARTS , NULL_CHARTS , NTEX";
    for (unsigned i = 0; i < textureStats.size(); ++i) {
        std::cout << " , " << "TEX" << i << "_SZ"
                  << " , " << "TEX" << i << "_INJ"
                  << " , " << "TEX" << i << "_OCCUPANCY";
    }
    std::cout << " , # INFO_TABLE_HEADER" << std::endl;

    std::cout << name << " , "
              << vn << " , " << fn << " , " << nonManifoldEdgeCount << " , " << nonManifoldVertexCount << " , "
              << atlasRegionCount << " , " << emptyRegionCount << " , " << textureNum;

    for (int i = 0; i < textureNum; ++i) {
        std::cout << " , " << textureStats[i].rw << "x" << textureStats[i].rh
                  << " , " << ((textureStats[i].overwrittenFragments == 0) ? "true" : "FALSE")
                  << " , " << (textureStats[i].totalFragments - textureStats[i].lostFragments) / ((double) textureStats[i].rw * textureStats[i].rh);
    }

    for (unsigned i = textureNum; i < textureStats.size(); ++i) {
        std::cout << " , , , ";
    }

    std::cout << " , # INFO_TABLE_ROW" << std::endl;

    std::cout << std::endl;

    LogAggregateStats("stats", graph, textureObject);

    ScaleTextureCoordinatesToImage(m, textureObject);

    double total_surface_area = 0;
    double total_uv_area = 0;
    for (auto& f : m.face) {
        double area3D = DistortionMetric::Area3D(f);
        double areaUV = std::abs(DistortionMetric::AreaUV(f));
        if (std::isfinite(area3D) && std::isfinite(areaUV)) {
            total_surface_area += area3D;
            total_uv_area += areaUV;
        }
    }


    for (auto& f : m.face) {
        // quality is the texel allocation per face area
        double area3D = DistortionMetric::Area3D(f);
        double areaUV = std::abs(DistortionMetric::AreaUV(f));
        if (std::isfinite(area3D) && std::isfinite(areaUV)) {
            double area = area3D * (total_uv_area / total_surface_area);
            double q = areaUV / area;
            ensure_condition(q >= 0);
            f.Q() = q;
        } else {
            f.Q() = 0;
        }

        if (tri::Index(m, f)%1000 == 0) {
            std::cout << area3D << " " << areaUV << " " << f.Q() << std::endl;
        }
    }

    std::cout << total_uv_area << " " << total_surface_area << std::endl;

    /*
    vcg::Histogramd h;
    tri::Stat<Mesh>::ComputePerFaceQualityHistogram(m, h);
    */

    vcg::Distribution<double> h;
    tri::Stat<Mesh>::ComputePerFaceQualityDistribution(m, h);

    double perc1 = h.Percentile(0.01);
    double perc99 = h.Percentile(0.99);

    std::cout << "Texel allocation per unit surface area range: " << perc1 << " - " << perc99 << std::endl;

    double rangeMin = 0;
    double rangeMax = 10;

    ensure_condition(perc1 >= rangeMin);
    ensure_condition(perc99 <= rangeMax);

    for (auto& f : m.face) {
        double q = f.Q();
        if (q < perc1)
            q = perc1;
        if (q > perc99)
            q = perc99;
        if (q < 1)
            f.C().lerp(vcg::Color4b(0, 100, 255, 255), vcg::Color4b(255, 48, 48, 255), 1.0 - q);
        else
            f.C().lerp(vcg::Color4b(0, 100, 255, 255), vcg::Color4b(48, 255, 48, 255), (q - 1) / (rangeMax - 1));

        //f.C().lerp(vcg::Color4b(0, 100, 255, 255), vcg::Color4b(255, 48, 48, 255), (q - perc1) / (perc99 - perc1));
        //f.C().lerp(vcg::Color4b(0, 100, 255, 255), vcg::Color4b(255, 48, 48, 255), (q - perc1) / (rangeMax - rangeMin));
    }

    m.textures.clear();
    tri::io::ExporterPLY<Mesh>::Save(m, "texel.ply", tri::io::Mask::IOM_FACEQUALITY | tri::io::Mask::IOM_FACECOLOR);

    GLTerminate();

    return 0;
}
