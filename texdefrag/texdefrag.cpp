#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/export.h>

#include <string>
#include <vector>
#include <iostream>

#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>

#include "mesh.h"
#include "uv.h"
#include "optimizer.h"
#include "texture_rendering.h"
#include "timer.h"

using namespace vcg;

int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " model minRegionSize(int) [--nofilter]" << std::endl;
        return -1;
    }

    bool filter = true;
    if (argc > 3 && std::string("--nofilter").compare(argv[3]) == 0) {
        std::cout << "Pull-Push filter disabled" << std::endl;
        filter = false;
    } else {
        std::cout << "Pull-Push filter enabled" << std::endl;
    }

    int minRegionSize = atoi(argv[2]);

    Mesh m;
    std::vector<std::shared_ptr<QImage>> imgVec;
    int loadMask;
    std::string modelName;

    if (LoadMesh(m, argv[1], imgVec, loadMask, modelName) == false) {
        std::cout << "Failed to open mesh" << std::endl;
        std::exit(-1);
    }

    assert(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    if (minRegionSize > m.FN()) {
        std::cout << "Error: minFaceCount > m.FN()" << std::endl;
        std::exit(-1);
    }

    StoreWedgeTexCoordAsAttribute(m);

    float uvMeshBorder;
    auto graph = ComputeParameterizationGraph(m, imgVec, &uvMeshBorder);

    // Print original info
    PrintParameterizationInfo(graph);

    Timer t;

#ifdef OLD_OPTIMIZER
    ReduceTextureFragmentation(m, *graph, minRegionSize);
#else
    ReduceTextureFragmentation(m, graph, minRegionSize);
#endif

    std::cout << "Processing took " << t.TimeElapsed() << " seconds" << std::endl;

    std::shared_ptr<QImage> img = RenderTexture(m, imgVec, filter);
    img->save(m.textures[0].c_str(), 0, 100);

    auto graph2 = ComputeParameterizationGraph(m, imgVec, &uvMeshBorder);
    // Print optimized info
    PrintParameterizationInfo(graph2);
    tri::io::Exporter<Mesh>::Save(m, modelName.c_str(), tri::io::Mask::IOM_WEDGTEXCOORD);

    return 0;
}


