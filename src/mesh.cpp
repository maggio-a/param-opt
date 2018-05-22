#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <iostream>
#include <string>

#include <QFileInfo>
#include <QDir>
#include <QString>

#include "mesh.h"
#include "gl_utils.h"

bool LoadMesh(Mesh &m, const char *fileName, TextureObjectHandle& textureObject, int &loadMask)
{
    m.Clear();
    textureObject = std::make_shared<TextureObject>();
    loadMask = 0;

    QFileInfo fi(fileName);

    if (!fi.exists() || !fi.isReadable()) {
        std::cout << "Unable to read " << fileName << std::endl;
        return false;
    }

    m.name = fi.fileName().toStdString();

    QString wd = QDir::currentPath();
    QDir::setCurrent(fi.absoluteDir().absolutePath());

    int r = tri::io::Importer<Mesh>::Open(m, fi.fileName().toStdString().c_str(), loadMask);
    if (r) {
        std::cout << tri::io::Importer<Mesh>::ErrorMsg(r) << std::endl;
        return false;
    }

    std::cout << "[LOG] Loaded mesh " << fileName << " (VN " <<  m.VN() << ", FN " << m.FN() << ")" << std::endl;

    for (const string& textureName : m.textures) {
        QFileInfo textureFile(textureName.c_str());
        textureFile.makeAbsolute();
        auto imgptr = std::make_shared<QImage>(textureFile.absoluteFilePath());
        if (!textureFile.exists() || !textureFile.isReadable() || imgptr->isNull()) {
            std::cout << "Unable to load texture file " << textureName.c_str() << std::endl;
            return false;
        }
        textureObject->AddImage(imgptr);
    }

    QDir::setCurrent(wd);
    return true;
}

bool SaveMesh(Mesh &m, const char *fileName, TextureObjectHandle& textureObject, bool color)
{
    int mask = tri::io::Mask::IOM_WEDGTEXCOORD;

    // Quick and dirty, make sure the texture extension is consistent
    for (std::string& textureName : m.textures) {
        textureName = textureName.substr(0, textureName.find_last_of('.')).append(".png");
    }

    if (color) mask = mask | tri::io::Mask::IOM_FACEQUALITY | tri::io::Mask::IOM_FACECOLOR;
    if (tri::io::Exporter<Mesh>::Save(m, fileName, mask)) {
        std::cout << "Error saving mesh file " << fileName << std::endl;
        return false;
    }

    QFileInfo fi(fileName);
    assert (fi.exists());

    QString wd = QDir::currentPath();
    QDir::setCurrent(fi.absoluteDir().absolutePath());

    for (std::size_t i = 0; i < textureObject->imgVec.size(); ++i) {
        if(textureObject->imgVec[i]->save(m.textures[i].c_str(), "png", 100) == false) {
            std::cout << "Error saving texture file " << m.textures[0] << std::endl;
            return false;
        }
    }

    QDir::setCurrent(wd);
    return true;
}

void BuildMeshFromFacePointers(Mesh &m, const std::vector<std::vector<Mesh::FacePointer>* >& vFpVecp)
{
    m.Clear();

    auto f = [&m](typename Mesh::FacePointer fptr) {
        tri::Allocator<Mesh>::AddFace(m, fptr->P(0), fptr->P(1), fptr->P(2));
    };

    for (auto fpVecp : vFpVecp) std::for_each(fpVecp->begin(), fpVecp->end(), f);

    tri::Clean<Mesh>::RemoveDuplicateVertex(m);
    tri::Allocator<Mesh>::CompactEveryVector(m);

    tri::UpdateTopology<Mesh>::FaceFace(m);
    tri::UpdateBounding<Mesh>::Box(m);
}

bool Parameterizable(Mesh &m)
{
    if (tri::Clean<Mesh>::CountNonManifoldEdgeFF(m) > 0) {
        return false;
    }

    if (tri::Clean<Mesh>::IsWaterTight(m)) {
        return false;
    }

    if (tri::Clean<Mesh>::MeshGenus(m) > 0) {
        return false;
    }

    return true;
}
