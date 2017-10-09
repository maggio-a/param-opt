#include "mesh.h"

#include <wrap/io_trimesh/import.h>

#include <iostream>
#include <string>

#include <QFileInfo>
#include <QDir>


bool LoadMesh(Mesh &m, const char *fileName, std::vector<std::shared_ptr<QImage>> &imgVec, int &loadMask, std::string &modelName)
{
    m.Clear();
    imgVec.clear();
    loadMask = 0;

    QFileInfo fi(fileName);

    if (!fi.exists() || !fi.isReadable()) {
        std::cout << "Unable to read " << fileName << std::endl;
        return false;
    }

    modelName = fi.fileName().toStdString();

    std::string wd = QDir::currentPath().toStdString();
    QDir::setCurrent(fi.absoluteDir().absolutePath());

    int r = tri::io::Importer<Mesh>::Open(m, fi.fileName().toStdString().c_str(), loadMask);
    if (r) {
        std::cout << tri::io::Importer<Mesh>::ErrorMsg(r) << std::endl;
        return false;
    }

    std::cout << "Loaded mesh " << fileName << " (VN " <<  m.VN() << ", FN " << m.FN() << ")" << std::endl;

    for (const string& textureName : m.textures) {
        auto imgptr = std::make_shared<QImage>(textureName.c_str());
        if (imgptr->isNull()) {
            std::cout << "Unable to load texture file " << textureName.c_str() << std::endl;
            return false;
        } else {
            std::cout << "Loaded texture " << textureName << std::endl;
        }
        imgVec.push_back(imgptr);
    }

    QDir::setCurrent(QString(wd.c_str()));
    return true;
}
