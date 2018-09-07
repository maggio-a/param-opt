#include "mesh.h"
#include "mesh_graph.h"
#include "uv.h"
#include "texture_optimization.h"
#include "texture_rendering.h"
#include "parameterization_checker.h"
#include "timer.h"
#include "gl_utils.h"
#include "texture.h"

#include <wrap/io_trimesh/io_mask.h>

#include <string>
#include <vector>
#include <iostream>

#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>



using namespace vcg;

int main(int argc, char *argv[])
{
    GLInit();

    Mesh m;
    TextureObjectHandle textureObject;
    int loadMask;

    assert(argc > 1);

    if (LoadMesh(m, argv[1], textureObject, loadMask) == false) {
        std::cout << "Failed to open mesh" << std::endl;
        std::exit(-1);
    }

    assert(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    tri::UpdateTopology<Mesh>::FaceFace(m);

    std::string name = m.name;
    int vn = m.VN();
    int fn = m.FN();

    int nonManifoldEdgeCount = tri::Clean<Mesh>::CountNonManifoldEdgeFF(m);
    int nonManifoldVertexCount = tri::Clean<Mesh>::CountNonManifoldVertexFF(m);

    int textureNum = textureObject->ArraySize();
    std::vector<std::pair<int,int>> textureSizes;
    for (int i = 0; i < textureNum; ++i)
        textureSizes.push_back(std::make_pair(textureObject->TextureWidth(i), textureObject->TextureHeight(i)));

    auto graph = ComputeParameterizationGraph(m, textureObject);

    std::size_t atlasRegionCount = graph->Count();

    std::size_t emptyRegionCount = 0;
    for (auto& entry : graph->charts)
        if (entry.second->AreaUV() == 0.0)
            emptyRegionCount++;

    std::vector<std::vector<Mesh::FacePointer>> facesByTexture(textureNum);
    for (auto& f : m.face) {
        int unit = f.WT(0).N();
        assert(unit < textureNum);
        facesByTexture[unit].push_back(&f);
    }

    std::vector<RasterizedParameterizationStats> textureStats(textureNum);
    for (int i = 0; i < textureNum; ++i) {
        textureStats[i] = GetRasterizationStats(facesByTexture[i], textureSizes[i].first, textureSizes[i].second);
    }

    std::cout << "NAME , VN , FN , NMANIF_EDGE , NMANIF_VERT , CHARTS , NULL_CHARTS , NTEX";
    for (int i = 0; i < textureNum; ++i) {
        std::cout << " , " << "TEX" << i << "_SZ"
                  << " , " << "TEX" << i << "_INJ"
                  << " , " << "TEX" << i << "_OCCUPANCY";
    }
    std::cout << " , # TABLE_INFO_HEADER" << std::endl;

    std::cout << name << " , "
              << vn << " , " << fn << " , " << nonManifoldEdgeCount << " , " << nonManifoldVertexCount << " , "
              << atlasRegionCount << " , " << emptyRegionCount << " , " << textureNum;

    for (int i = 0; i < textureNum; ++i) {
        std::cout << " , " << textureSizes[i].first << "x" << textureSizes[i].second
                  << " , " << ((textureStats[i].overwrittenFragments > 0) ? "true" : "FALSE")
                  << " , " << (textureStats[i].totalFragments - textureStats[i].overwrittenFragments) / ((double) textureStats[i].rw * textureStats[i].rh);
    }

    std::cout << " , # TABLE_INFO_ENTRY" << std::endl;

    std::cout << std::endl;


    GLTerminate();

    return 0;
}
