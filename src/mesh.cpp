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
#include "timer.h"
#include "utils.h"

bool LoadMesh(Mesh &m, const char *fileName, TextureObjectHandle& textureObject, int &loadMask)
{
    m.Clear();
    textureObject = std::make_shared<TextureObject>();
    loadMask = 0;

    QFileInfo fi(fileName);
    fi.makeAbsolute();

    if (!fi.exists() || !fi.isReadable()) {
        std::cout << "Unable to read " << fileName << std::endl;
        return false;
    }

    std::string dirname = fi.dir().dirName().toStdString();
    m.name = dirname + "_" + fi.fileName().toStdString();

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
        if (!textureFile.exists() || !textureFile.isReadable()) {
            std::cout << "Error: texture file " << textureName.c_str() << " does not exist or is not readable." << std::endl;
            return false;
        }

        auto imgptr = std::make_shared<QImage>(textureFile.absoluteFilePath());
        if (imgptr->isNull()) {
            std::cout << "Error: failed to load texture file " << textureName.c_str() << std::endl;
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
    for (std::size_t i = 0; i < m.textures.size(); ++i) {
        std::stringstream suffix;
        suffix << "_texture_" << i << ".png";
        std::string s(fileName);
        m.textures[i] = s.substr(0, s.find_last_of('.')).append(suffix.str());
    }

    std::cout << "Saving mesh file... ";
    if (color) mask = mask | tri::io::Mask::IOM_FACEQUALITY | tri::io::Mask::IOM_FACECOLOR;
    int err;
    if ((err = tri::io::Exporter<Mesh>::Save(m, fileName, mask))) {
        std::cout << "Error saving mesh file " << fileName << std::endl;
        std::cout << tri::io::Exporter<Mesh>::ErrorMsg(err) << std::endl;
        return false;
    }
    std::cout << " done." << std::endl;

    QFileInfo fi(fileName);
    ensure_condition (fi.exists());

    QString wd = QDir::currentPath();
    QDir::setCurrent(fi.absoluteDir().absolutePath());

    std::cout << "Saving texture files... ";
    Timer t;
    for (std::size_t i = 0; i < textureObject->imgVec.size(); ++i) {
        if (textureObject->imgVec[i]->save(m.textures[i].c_str(), "png", 66) == false) {
            std::cout << "Error saving texture file " << m.textures[i] << std::endl;
            return false;
        }
    }
    std::cout << " done." << std::endl;
    std::cout << "Writing textures took " << t.TimeElapsed() << " seconds" << std::endl;

    QDir::setCurrent(wd);
    return true;
}

void MeshFromFacePointers(const std::vector<Mesh::FacePointer>& vfp, Mesh& out)
{
    out.Clear();
    std::unordered_map<Mesh::VertexPointer, Mesh::VertexPointer> vpmap;
    vpmap.reserve(vfp.size() * 2);
    std::size_t vn = 0;
    for (auto fptr : vfp) {
        for (int i = 0; i < 3; ++i) {
            if (vpmap.count(fptr->V(i)) == 0) {
                vn++;
                vpmap[fptr->V(i)] = nullptr;
            }
        }
    }
    auto mvi = tri::Allocator<Mesh>::AddVertices(out, vn);
    auto mfi = tri::Allocator<Mesh>::AddFaces(out, vfp.size());
    for (auto fptr : vfp) {
        Mesh::FacePointer mfp = &*mfi++;
        for (int i = 0; i < 3; ++i) {
            Mesh::VertexPointer vp = fptr->V(i);
            typename Mesh::VertexPointer& mvp = vpmap[vp];
            if (mvp == nullptr) {
                mvp = &*mvi++;
                mvp->P() = vp->P();
            }
            mfp->V(i) = mvp;
            mfp->WT(i) = fptr->WT(i);
        }
        mfp->SetMesh();
    }
}

bool Parameterizable(Mesh &m)
{
    tri::UpdateTopology<Mesh>::FaceFace(m);

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
