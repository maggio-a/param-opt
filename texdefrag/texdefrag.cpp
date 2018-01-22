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
#include "mesh_graph.h"
#include "uv.h"
#include "texture_optimization.h"
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
    TextureObjectHandle textureObject;
    int loadMask;
    std::string modelName;

    if (LoadMesh(m, argv[1], textureObject, loadMask, modelName) == false) {
        std::cout << "Failed to open mesh" << std::endl;
        std::exit(-1);
    }

    assert(textureObject->ArraySize() == 1 && "Currently only single texture is supported");

    assert(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    if (minRegionSize > m.FN()) {
        std::cout << "WARNING: minFaceCount > m.FN()" << std::endl;
    }

    StoreWedgeTexCoordAsAttribute(m);

    float uvMeshBorder;
    auto graph = ComputeParameterizationGraph(m, textureObject, &uvMeshBorder);

    // Print original info
    PrintParameterizationInfo(graph);

    Timer t;

    GraphManager gm{graph};

    ReduceTextureFragmentation_NoPacking(gm, minRegionSize);
/*
    std::cout << "Rendering texture..." << std::endl;
    TextureObjectHandle newTexture = RenderTexture(m, textureObject, filter, nullptr);

    std::cout << "Processing took " << t.TimeElapsed() << " seconds" << std::endl;

    if (SaveMesh(m, modelName.c_str(), newTexture) == false) {
        std::cout << "Model not saved correctly" << std::endl;
    }

    auto graph2 = ComputeParameterizationGraph(m, textureObject, &uvMeshBorder);
    // Print optimized info
    PrintParameterizationInfo(graph2);
*/
    return 0;
}


