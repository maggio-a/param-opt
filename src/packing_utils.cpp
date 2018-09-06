#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/outline_support.h>
#include <wrap/qt/outline2_rasterizer.h>

#include <vector>

#include "packing_utils.h"

using namespace vcg;

/* Returns a scaled version of the outline such that it lies in [0,1]x[0,1] */
Outline2f NormalizeOutline(const Outline2f& outline)
{
    Outline2f normalized = outline;
    Box2f bb;
    for (const auto& p : outline)
        bb.Add(p);
    float scale = std::max(bb.DimX(), bb.DimY());
    for (auto& p : normalized)
        p /= scale;
    return normalized;
}

/* Returns the normalized difference between the area computed as the top-bottom
 * horizons and the actual rasterization area, in absolute value. This measure is
 * averaged over 8 rotations of the outline */
float EvaluatePackingQuality_Horizon(const Outline2f& outline)
{
    constexpr int rotationNum = 8;
    constexpr float scale = 128.0f;
    Outline2f normalizedOutline = NormalizeOutline(outline);
    RasterizedOutline2 rasterization;
    rasterization.resetState(rotationNum);
    rasterization.setPoints(normalizedOutline);
    for (int i = 0; i < rotationNum / 4; ++i)
        QtOutline2Rasterizer::rasterize(rasterization, scale, i, rotationNum, 1, 0);

    float qualityIndex = 0;
    for (int k = 0 ; k < rotationNum; ++k) {
        const std::vector<std::vector<int>>& grid = rasterization.getGrids(k);
        const std::vector<int>& left = rasterization.getLeft(k);
        const std::vector<int>& deltaX = rasterization.getDeltaX(k);
        int area = 0;
        for (unsigned i = 0; i < left.size(); ++i) {
            for (int j = left[i]; j < left[i] + deltaX[i]; ++j)
                if (grid[i][j] > 0)
                    area++;
        }
        double gridArea = rasterization.gridWidth(k) * rasterization.gridHeight(k);
        qualityIndex += (rasterization.getDiscreteArea(k) - area) / area;
    }

    return 1.0f - (qualityIndex / rotationNum);
}

float EvaluatePackingQuality_VoidArea(const Outline2f& outline)
{
    /*
    constexpr int rotationNum = 8;
    constexpr float scale = 128.0f;
    Outline2f normalizedOutline = NormalizeOutline(outline);
    RasterizedOutline2 rasterization;
    rasterization.resetState(rotationNum);
    rasterization.setPoints(normalizedOutline);
    for (int i = 0; i < rotationNum / 4; ++i)
        QtOutline2Rasterizer::rasterize(rasterization, scale, i, rotationNum, 1, 1);

    float qualityIndex = 0;
    for (int k = 0 ; k < rotationNum; ++k) {
        const std::vector<std::vector<int>>& grid = rasterization.getGrids(k);
        const std::vector<int>& left = rasterization.getLeft(k);
        const std::vector<int>& deltaX = rasterization.getDeltaX(k);
        int area = 0;
        for (unsigned i = 0; i < left.size(); ++i) {
            for (int j = left[i]; j < left[i] + deltaX[i]; ++j)
                if (grid[i][j] > 0)
                    area++;
        }
        double gridArea = rasterization.gridWidth(k) * rasterization.gridHeight(k);
        qualityIndex += (gridArea - area) / gridArea;
    }

    return (qualityIndex / rotationNum);
    */
    constexpr int rotationNum = 100;
    Outline2f normalizedOutline = NormalizeOutline(outline);

    float qualityIndex = std::numeric_limits<float>::max();
    float rot = (2.0 * M_PI) / rotationNum;
    for (int k = 0 ; k < rotationNum; ++k) {
        vcg::Box2f bb;
        for (auto& p : normalizedOutline) {
            p.Rotate(rot);
            bb.Add(p);
        }
        float outlineArea = std::abs(tri::OutlineUtil<float>::Outline2Area(normalizedOutline));
        float boxArea = bb.DimX() * bb.DimY();
        float voidArea = boxArea - outlineArea;
        qualityIndex = std::min(voidArea / boxArea, qualityIndex);
    }

    return 1.0f - qualityIndex; // returns the best filled fraction
}

float EvaluatePackingQuality_Perimeter(const Outline2f& outline)
{
    constexpr int rotationNum = 100;
    Outline2f normalizedOutline = NormalizeOutline(outline);

    float qualityIndex = std::numeric_limits<float>::max();
    float rot = (2.0 * M_PI) / rotationNum;
    for (int k = 0 ; k < rotationNum; ++k) {
        vcg::Box2f bb;
        for (auto& p : normalizedOutline) {
            p.Rotate(rot);
            bb.Add(p);
        }
        float perimeter = tri::OutlineUtil<float>::Outline2Perimeter(normalizedOutline);
        float boxPerimeter = 2 * bb.DimX() + 2 * bb.DimY();
        qualityIndex = std::min(std::abs(perimeter - boxPerimeter), qualityIndex);
    }

    return 1.0f / (qualityIndex + 1e-8);
}

