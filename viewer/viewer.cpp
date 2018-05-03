#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/io_mask.h>

#include <wrap/io_trimesh/export.h>

#include <string>
#include <vector>
#include <iostream>

#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>

#include "mesh.h"
#include "mesh_viewer.h"
#include "texture_optimization.h"
#include "uv.h"
#include "gl_utils.h"


using namespace vcg;

int main(int argc, char *argv[])
{
     if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " model minRegionSize(int) [--nofilter]" << std::endl;
        std::exit(-1);
    }

    int minRegionSize = atoi(argv[2]);
    assert(minRegionSize > 0);

    Mesh m;
    TextureObjectHandle textureObject;
    int loadMask;
    std::string fileName;

    if (LoadMesh(m, argv[1], textureObject, loadMask, fileName) == false) {
        std::cout << "Failed to open mesh" << std::endl;
        std::exit(-1);
    }

    //assert(textureObject->ArraySize() == 1 && "Currently only single texture is supported");

    assert(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    ComputeParameterizationScaleInfo(m);

    StoreWedgeTexCoordAsAttribute(m);
    PreprocessMesh(m);

    float uvMeshBorder;
    auto graph = ComputeParameterizationGraph(m, textureObject, &uvMeshBorder);

    // Print original info
    PrintParameterizationInfo(graph);

    GLInit();

    MeshViewer viewer(graph, std::size_t(minRegionSize), fileName);
    viewer.Run();

    GLTerminate();

    return 0;
}
