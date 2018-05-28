#include <wrap/io_trimesh/io_mask.h>

#include <string>
#include <vector>
#include <iostream>

#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>

#include "mesh.h"
#include "mesh_attribute.h"
#include "texture_rendering.h"
#include "gl_utils.h"

using namespace vcg;

int main(int argc, char *argv[])
{
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " model_src model_dest" << std::endl;
        std::exit(-1);
    }

    Mesh m1;
    TextureObjectHandle t1;
    int loadMask;

    Mesh m2;
    TextureObjectHandle t2;

    if (LoadMesh(m1, argv[1], t1, loadMask) == false) {
        std::cout << "Failed to open mesh " << argv[1] << std::endl;
        std::exit(-1);
    }
    assert(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);
    if (LoadMesh(m2, argv[2], t2, loadMask) == false) {
        std::cout << "Failed to open mesh " << argv[2] << std::endl;
        std::exit(-1);
    }
    assert(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    auto wtcsattr = GetWedgeTexCoordStorageAttribute(m2);

    assert(m1.FN() == m2.FN());

    t2 = nullptr; // deallocate unused image data

    for (int i = 0; i < m1.FN(); ++i) {
        auto& f1 = m1.face[i];
        auto& f2 = m2.face[i];

        for (int j = 0; j < 3; ++j) {
            wtcsattr[f2].tc[j] = f1.WT(j);
        }
    }

    m1.Clear();

    GLInit();

    TextureObjectHandle newTexture = RenderTexture(m2, t1, true, InterpolationMode::Linear, nullptr);
    std::string savename = "out_" + m1.name;
    if (SaveMesh(m2, savename.c_str(), newTexture, true) == false) {
        std::cout << "Model not saved correctly" << std::endl;
    }

    GLTerminate();

    return 0;
}
