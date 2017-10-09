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
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " model" << std::endl;
        std::exit(-1);
    }

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

    MeshViewer viewer(graph);
    viewer.Run();

    return 0;

}
