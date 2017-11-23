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
    TextureObjectHandle textureObject;
    int loadMask;
    std::string fileName;

    if (LoadMesh(m, argv[1], textureObject, loadMask, fileName) == false) {
        std::cout << "Failed to open mesh" << std::endl;
        std::exit(-1);
    }

    //assert(textureObject->ArraySize() == 1 && "Currently only single texture is supported");

    assert(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    StoreWedgeTexCoordAsAttribute(m);

    float uvMeshBorder;
    auto graph = ComputeParameterizationGraph(m, textureObject, &uvMeshBorder);

    // Print original info
    PrintParameterizationInfo(graph);

    MeshViewer viewer(graph, std::size_t(minRegionSize), fileName);
    viewer.Run();

    return 0;

    for (auto &f : m.face) {
        for (int i =0 ; i < 3; ++i) {
            f.V(i)->T() = f.WT(i);
        }
    }

/*
    tri::UpdateTopology<Mesh>::FaceFace(m);
    tri::UpdateFlags<Mesh>::VertexBorderFromFaceAdj(m);
    tri::AreaPreservingTextureOptimizer<Mesh> opt(m);

    opt.TargetCurrentGeometry();
    opt.SetBorderAsFixed();

    for (int i = 0; i < 100; ++i) {
        Timer t;
        opt.Iterate();
        tri::io::ExporterOBJ<Mesh>::Save(m, "opt.obj", tri::io::Mask::IOM_WEDGTEXCOORD);
        tri::UpdateTexture<Mesh>::WedgeTexFromVertexTex(m);
        std::cout << "Iteration took " << t.TimeSinceLastCheck() << " seconds";
        std::cout << std::endl;
    }

    return 0; */

}
