#ifndef PACKING_UTILS_H
#define PACKING_UTILS_H

#include <vcg/complex/algorithms/outline_support.h>
#include <vcg/space/rasterized_outline2_packer.h>

#include <fstream>

#include "utils.h"
#include "logging.h"


class Mesh;

using Outline2f = std::vector<vcg::Point2f>;

enum PackingScore {
    HorizonBased,
    VoidBased,
    PerimeterAreaRatio
};

float ComputeAtlasScore(const std::vector<Outline2f>& outlines, PackingScore ps);

float EvaluatePackingQuality_Horizon(const Outline2f& outline);
float EvaluatePackingQuality_VoidArea(const Outline2f& outline);
float EvaluatePackingQuality_Perimeter(const Outline2f& outline);


template <class ScalarType, class RasterizerType>
class PixelShiftingOptimizer
{
    using RasterizedOutline2 = vcg::RasterizedOutline2;

    typedef typename vcg::Box2<ScalarType> Box2x;
    typedef typename vcg::Point2<ScalarType> Point2x;
    typedef typename vcg::Similarity2<ScalarType> Similarity2x;
    typedef typename std::vector<Point2x> Poly;

    static inline constexpr unsigned PIXEL_CLEAR()
    {
        return 0xFFFFFFFF;
    }

    static inline bool COLLIDE(const std::vector<std::vector<unsigned>>& buffer, int row, int col, unsigned k)
    {
        return (buffer[row][col] != PIXEL_CLEAR()) && (buffer[row][col] != k);
    }

    static bool Collision(const std::vector<std::vector<unsigned>>& buffer,
                          const std::vector<std::vector<vcg::Point2i>>& boundaries,
                          int xpos, int ypos,
                          unsigned k)
    {
        for (const auto& p : boundaries[k]) {
            int row = ypos + p.Y();
            int col = xpos + p.X();
            if (COLLIDE(buffer, row, col, k))
                return true;
        }
        return false;
    }

    static void Set(std::vector<std::vector<unsigned>>& buffer,
                    const std::vector<std::vector<vcg::Point2i>>& boundaries,
                    int xpos, int ypos,
                    unsigned k)
    {
        for (const auto& p : boundaries[k]) {
            int row = ypos + p.Y();
            int col = xpos + p.X();
            //ensure_condition(!COLLIDE(buffer, row, col, k));
            bool collision = COLLIDE(buffer, row, col, k);
            if (collision) {
                LOG_DEBUG << "collision at " << col << " " << row;
                //ensure_condition(!collision);
            }

            buffer[row][col] = k;
        }
    }

    static void Clear(std::vector<std::vector<unsigned>>& buffer,
                      const std::vector<std::vector<vcg::Point2i>>& boundaries,
                      int xpos, int ypos,
                      unsigned k)
    {
        for (const auto& p : boundaries[k]) {
            int row = ypos + p.Y();
            int col = xpos + p.X();
            //ensure_condition(!COLLIDE(buffer, row, col, k));
            bool collision = COLLIDE(buffer, row, col, k);
            if (collision) {
                LOG_DEBUG << "collision at " << col << " " << row;
                ensure_condition(!collision);
            }
            buffer[row][col] = PIXEL_CLEAR();
        }
    }

    static bool IsBoundary(const std::vector<std::vector<int>>& grid, int row, int col, int nr, int nc)
    {
        if (row == 0 || row == nr - 1)
            return true;
        if (col == 0 || col == nc - 1)
            return true;
        for (int i = -1; i < 2; ++i) {
            for (int j = -1; j < 2; ++j) {
                if (grid[row + i][col + j] == 0)
                    return true;
            }
        }
        return false;
    }

    static void ComputeBoundaryFromGrid(const std::vector<std::vector<int>>& grid, std::vector<vcg::Point2i>& boundaryVector)
    {
        ensure_condition(grid.size() > 0);
        int numRows = grid.size();
        ensure_condition(grid[0].size() > 0);
        int numCols = grid[0].size();
        boundaryVector.clear();
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                if (grid[i][j] == 1 && IsBoundary(grid, i, j, numRows, numCols)) {
                    boundaryVector.push_back(vcg::Point2i(j, i));
                }
            }
        }
    }

public:

    struct Parameters {
        int w;
        int h;
        int gutter;
    };

    static vcg::Point2i Apply(std::vector<Poly>& polyVec, std::vector<Similarity2x>& trVec, const Parameters& param)
    {
        LOG_INFO << "Applying pixel shifting";
        using vcg::Point2i;

        // compute scale factor that optimally fits the raster space
        ScalarType scale = std::min(ScalarType(param.w), ScalarType(param.h));
        ensure_condition(scale > 0);

        LOG_VERBOSE << "Grid size is " << param.w << "x" << param.h;

        // compute the initial offsets
        std::vector<Box2x> polyBoxes;
        for (const auto& poly : polyVec) {
            Box2x polyBox;
            for (const auto& point : poly)
                polyBox.Add(point);
            polyBoxes.push_back(polyBox);
        }


        LOG_VERBOSE << "Rasterizing outlines... ";

        // initialize the boundary buffers
        std::vector<std::vector<Point2i>> boundaries(polyVec.size());
        std::vector<Point2i> extents(polyVec.size());
        for (unsigned k = 0; k < polyVec.size(); ++k) {
            RasterizedOutline2 rasterization;
            rasterization.resetState(4);
            rasterization.setPoints(polyVec[k]);
            RasterizerType::rasterize(rasterization, scale, 0, 4, 1, param.gutter); // needs 4 rot for compatibility with the outline2_rasterizer
            ComputeBoundaryFromGrid(rasterization.getGrids(0), boundaries[k]);
            extents[k] = Point2i(rasterization.gridWidth(0), rasterization.gridHeight(0));
            /*
            std::vector<std::vector<int>> grid;
            RasterizerType::rasterize(polyVec[k], scale, 0, param.gutter, grid);
            ComputeBoundaryFromGrid(grid, boundaries[k]);
            extents[k] = Point2i(grid[0].size(), grid.size());
            */
        }

        // initialize the raster buffer with the object outlines
        const int margin = 32;
        std::vector<std::vector<unsigned>> buffer(param.h + margin);
        for (auto& row : buffer)
            row.resize(param.w + margin, PIXEL_CLEAR());

        std::vector<int> xpos(polyVec.size(), -1);
        std::vector<int> ypos(polyVec.size(), -1);

        LOG_VERBOSE << "Initializing buffer...";
        for (unsigned k = 0; k < polyVec.size(); ++k) {
            xpos[k] = std::max(int(std::round(polyBoxes[k].min.X() * scale)), 0);
            ypos[k] = std::max(int(std::round(polyBoxes[k].min.Y() * scale)), 0);
            //xpos[k] = std::max(int(std::round(polyBoxes[k].min.X() * scale)) - param.gutter, 0);
            //ypos[k] = std::max(int(std::round(polyBoxes[k].min.Y() * scale)) - param.gutter, 0);

            vcg::Box2i extent;
            for (auto& p : boundaries[k])
                extent.Add(p);

            Set(buffer, boundaries, xpos[k], ypos[k], k);
        }

        /*
        {
            std::ofstream ofs("shift.ppm", std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
            if (!ofs) {
                ensure_condition(0);
            } else {
                ofs << "P6 " << param.w << " " << param.h << " 255\n";
                for (int i = 0; i < param.h; ++i) {
                    for (int j = 0; j < param.w; ++j) {
                        if (buffer[i][j] == PIXEL_CLEAR()) {
                            ofs << (unsigned char) 255 << (unsigned char) 255 << (unsigned char) 255;
                        } else {
                            ofs << (unsigned char) 255 << (unsigned char) 0 << (unsigned char) 0;
                        }
                    }
                }
                ofs.close();
            }
        }
        */

        LOG_VERBOSE << "Shifting charts...";

        bool bufferChanged = true;
        int count = 0;
        while (bufferChanged) {
            bufferChanged = false;
            for (unsigned k = 0; k < polyVec.size(); ++k) {
                int y = ypos[k];
                int x = xpos[k];
                if (x > 0 && y > 0 && !Collision(buffer, boundaries, x - 1, y - 1, k)) {
                    x--;
                    y--;
                } else if (y > 0 && !Collision(buffer, boundaries, x, y - 1, k)) {
                    y--;
                } else if (x > 0 && !Collision(buffer, boundaries, x - 1, y, k)) {
                    x--;
                }
                if (x != xpos[k] || y != ypos[k]) {
                    Clear(buffer, boundaries, xpos[k], ypos[k], k);
                    Set(buffer, boundaries, x, y, k);
                    xpos[k] = x;
                    ypos[k] = y;
                    bufferChanged = true;
                }
            }
            count++;
        }

        vcg::Box2i cover;
        for (unsigned k = 0; k < polyVec.size(); ++k) {
            for (const auto& p : boundaries[k])
                cover.Add(Point2i(xpos[k], ypos[k]) + p);
        }

        // compute the vector of transforms
        trVec.clear();
        for (unsigned k = 0; k < polyVec.size(); ++k) {
            Similarity2x sim;

            ScalarType offsetX = (extents[k].X() - ceil(scale * polyBoxes[k].DimX())) / 2.0;
            ScalarType offsetY = (extents[k].Y() - ceil(scale * polyBoxes[k].DimY())) / 2.0;

            //ScalarType offsetX = 1;
            //ScalarType offsetY = 1;

            sim.rotRad = 0;
            sim.sca = scale;
            //Point2x finalOffset(ScalarType(xpos[k]), ScalarType(ypos[k]));
            //sim.tra = - (scale * initialOffset[k]) + finalOffset;
            sim.tra = vcg::Point2<ScalarType>(xpos[k] - scale*polyBoxes[k].min.X() + offsetX,
                                         ypos[k] - scale*polyBoxes[k].min.Y() + offsetY);
                                         //param.h - (ypos[k] + extents[k].Y()) - scale*polyBoxes[k].min.Y() + offsetY);
            trVec.push_back(sim);
        }

        /*
        std::ofstream ofs("shift_final.ppm", std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
        if (!ofs) {
            ensure_condition(0);
        } else {
            ofs << "P6 " << param.w << " " << param.h << " 255\n";
            for (int i = 0; i < param.h; ++i) {
                for (int j = 0; j < param.w; ++j) {
                    if (buffer[i][j] == PIXEL_CLEAR()) {
                        ofs << (unsigned char) 255 << (unsigned char) 255 << (unsigned char) 255;
                    } else {
                        ofs << (unsigned char) 255 << (unsigned char) 0 << (unsigned char) 0;
                    }
                }
            }
            ofs.close();
        }
        */

        LOG_INFO << "Pixel shifting finished";

        double aspectRatio;
        if (param.w == param.h) {
            int maxDim = std::max(cover.DimX(), cover.DimY());
            return Point2i(maxDim, maxDim);
        } else {
            aspectRatio = param.w / (double) param.h;
            if (aspectRatio > 1) {
                ensure_condition(cover.DimX() >= cover.DimY());
                return Point2i(cover.DimX(), cover.DimX() / aspectRatio);
            } else {
                ensure_condition(cover.DimY() >= cover.DimX());
                return Point2i(cover.DimY() * aspectRatio, cover.DimY());
            }
        }
    }

};


#endif // PACKING_UTILS_H

