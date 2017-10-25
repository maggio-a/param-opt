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
#include "mesh_viewer.h"
#include "uv.h"

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
    std::vector<std::shared_ptr<QImage>> imgVec;
    int loadMask;
    std::string modelName;

    if (LoadMesh(m, argv[1], imgVec, loadMask, modelName) == false) {
        std::cout << "Failed to open mesh" << std::endl;
        std::exit(-1);
    }

    assert(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    StoreWedgeTexCoordAsAttribute(m);

    float uvMeshBorder;
    auto graph = ComputeParameterizationGraph(m, imgVec, &uvMeshBorder);

    // Print original info
    PrintParameterizationInfo(graph);

    MeshViewer viewer(graph, std::size_t(minRegionSize));
    viewer.Run();

    return 0;

}
