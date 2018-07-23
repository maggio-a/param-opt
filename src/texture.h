#ifndef TEXTURE_H
#define TEXTURE_H

#include "mesh.h"

#include "mesh_graph.h"
#include "gl_utils.h"

#include <vector>

#include <QImage>
#include <QPainter>
#include <vcg/space/rect_packer.h>

void GenerateDistortionTextures(Mesh& m, TextureObjectHandle textureObject);

void EvaluateGutterResistance(Mesh& m, TextureObjectHandle textureObject);

static void CompactTextureData(Mesh& m, TextureObjectHandle texture, double *coverage)
{
    if (texture->ArraySize() == 1) {
        *coverage = 1.0;
        return;
    }

    int ntex = texture->imgVec.size();
    const int minRes = 1024;

    std::vector<Point2i> sizes;
    long texels = 0;
    for (auto qimgptr : texture->imgVec) {
        int w = qimgptr->width();
        int h = qimgptr->height();
        assert((w % minRes) == 0 && (h % minRes) == 0);
        sizes.push_back(Point2i{w / minRes, h / minRes});
        texels += w * h;
    }

    Point2i grid{16, 16};
    std::vector<Point2i> positions;
    Point2i cover;
    bool packed = RectPacker<double>::PackInt(sizes, grid, positions, cover);
    assert(packed);

    tri::Allocator<Mesh>::DeletePerMeshAttribute(m, "materialVector");
    tri::Allocator<Mesh>::DeletePerFaceAttribute(m, "materialIndex");

    std::string texName = std::string("out_") + m.textures[0];
    m.textures.clear();
    m.textures.push_back(texName);

    std::vector<Point2d> scale;
    std::vector<Point2d> translate;
    for (int i = 0; i < ntex; ++i) {
        double scale_u = sizes[i][0] / (double) grid[0];
        double scale_v = sizes[i][1] / (double) grid[1];
        double u0 = positions[i][0] / (double) grid[0];
        double v0 = 1.0 - ((positions[i][1] + sizes[i][1]) / (double) grid[1]);
        scale.push_back(Point2d{scale_u, scale_v});
        translate.push_back(Point2d{u0, v0});
    }

    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            int n = f.WT(i).N();
            f.WT(i).P().Scale(scale[n][0], scale[n][1]);
            f.WT(i).P() += translate[n];
            f.WT(i).N() = 0;
        }
    }

    std::shared_ptr<QImage> img = std::make_shared<QImage>(16384, 16384, texture->imgVec[0]->format());
    QPainter painter{img.get()};
    for (int i = 0; i < ntex; ++i) {
        Point2i p = minRes * positions[i];
        const QImage& imread = *(texture->imgVec[i]);
        painter.drawImage(p[0], p[1], imread);
    }
    painter.end();

    assert(img->save("test.jpg", 0, 100));
    //tri::io::Exporter<Mesh>::Save(m, "test.obj", tri::io::Mask::IOM_WEDGTEXCOORD);

    *coverage = texels / (double) (img->width() * img->height());
    texture->imgVec.clear();
    texture->AddImage(img);
}

#endif // TEXTURE_H

