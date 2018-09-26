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

    LogAggregateStats(m.name, graph, textureObject);

    GLTerminate();

    return 0;
}
