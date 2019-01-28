#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/outline_support.h>
#include <wrap/qt/outline2_rasterizer.h>

#include <vector>

#include "packing_utils.h"


/* Returns a scaled version of the outline such that its area is unitary */
static Outline2f NormalizeOutline(const Outline2f& outline);


using namespace vcg;

std::vector<Outline2f> NormalizeAtlas(const std::vector<Outline2f>& outlines)
{
    float areasum = 0;
    for (auto& outline : outlines) {
        areasum += tri::OutlineUtil<float>::Outline2Area(outline);
    }

    float scale = 1.0f / std::sqrt(areasum);
    std::vector<Outline2f> normalizedAtlas = outlines;
    for (auto& o : normalizedAtlas)
        for (auto& po : o)
            po *= scale;

    return normalizedAtlas;
}

float ComputeAtlasScore(const std::vector<Outline2f>& outlines, PackingScore ps)
{
    //std::vector<Outline2f> outlines = NormalizeAtlas(outlines_);

    float scoresum = 0;
    float areasum = 0;

    for (auto& outline : outlines) {
        float area = tri::OutlineUtil<float>::Outline2Area(outline);
        ensure_condition(std::isfinite(area));
        float score;
        switch (ps) {
        case PackingScore::HorizonBased:
            score = EvaluatePackingQuality_Horizon(outline);
            break;
        case PackingScore::VoidBased:
            score = EvaluatePackingQuality_VoidArea(outline);
            break;
        case PackingScore::PerimeterAreaRatio:
            score = EvaluatePackingQuality_Perimeter(outline);
            break;
        default:
            score = -1;
        }
        scoresum += area * score;
        areasum += area;
    }

    return (scoresum / areasum);
}

static Outline2f NormalizeOutline(const Outline2f& outline)
{
    return outline;
    /*
    Outline2f normalized = outline;
    float area = tri::OutlineUtil<float>::Outline2Area(outline);
    float sf = std::sqrt(area);
    Box2f bb;
    for (auto& p : normalized) {
        p /= sf;
        bb.Add(p);
    }
    float newarea = tri::OutlineUtil<float>::Outline2Area(normalized);
    LOG_INFO << "Outline normalization : " << newarea << " | bb = " << bb.DimX() << " " << bb.DimY();
    return normalized;
    */
}

/* Returns the normalized difference between the area between the top and bottom
 * horizons and the actual rasterization area, in absolute value. This measure is
 * averaged over 8 rotations of the outline.
 * The minimum of this index is one and is attained when the shapes are compact
 * (e.g. circles and squares).
 */
float EvaluatePackingQuality_Horizon(const Outline2f& outline)
{
    constexpr int rotationNum = 8;
    constexpr float scale = 1.0f;
    Outline2f normalizedOutline = NormalizeOutline(outline);
    RasterizedOutline2 rasterization;
    rasterization.resetState(rotationNum);
    rasterization.setPoints(normalizedOutline);
    for (int i = 0; i < rotationNum / 4; ++i)
        QtOutline2Rasterizer::rasterize(rasterization, scale, i, rotationNum, 1, 0);

    double qualityIndex = 0;
    for (int k = 0 ; k < rotationNum; ++k) {
        const std::vector<std::vector<int>>& grid = rasterization.getGrids(k);
        //const std::vector<int>& left = rasterization.getLeft(k);
        //const std::vector<int>& deltaX = rasterization.getDeltaX(k);
        int area = 0;

        for (unsigned i = 0; i < grid.size(); ++i) {
            for (unsigned j = 0; j < grid[i].size(); ++j)
                if (grid[i][j] > 0)
                    area++;
        }

        /*
        for (unsigned i = 0; i < left.size(); ++i) {
            for (int j = left[i]; j < left[i] + deltaX[i]; ++j)
                if (grid[i][j] > 0)
                    area++;
        }
        */

        double gridArea = rasterization.gridWidth(k) * rasterization.gridHeight(k);
        if (rasterization.getDiscreteArea(k) < area) {
            LOG_ERR << rasterization.getDiscreteArea(k) << "    " << area;
        }
        qualityIndex += (rasterization.getDiscreteArea(k) - area) / (double) area;
    }

    return (float) (qualityIndex / rotationNum);
}

/* Returns the normalized difference between the bounding box area and the
 * outline area. The maximum is 1 and is attained when the shape is compact. */
float EvaluatePackingQuality_VoidArea(const Outline2f& outline)
{
    constexpr int rotationNum = 8;
    Outline2f normalizedOutline = NormalizeOutline(outline);
    float outlineArea = std::abs(tri::OutlineUtil<float>::Outline2Area(normalizedOutline));
    float qualityIndex = 0;
    float rot = (2.0 * M_PI) / rotationNum;
    for (int k = 0 ; k < rotationNum; ++k) {
        vcg::Box2f bb;
        for (auto& p : normalizedOutline) {
            p.Rotate(rot);
            bb.Add(p);
        }
        float boxArea = bb.DimX() * bb.DimY();
        float voidArea = boxArea - outlineArea;
        qualityIndex += (voidArea / boxArea);
    }

    return 1.0f - (qualityIndex / rotationNum);
}

/* Returns the squared perimeter to area ratio */
float EvaluatePackingQuality_Perimeter(const Outline2f& outline)
{
    Outline2f normalizedOutline = NormalizeOutline(outline);
    float perimeter = tri::OutlineUtil<float>::Outline2Perimeter(normalizedOutline);
    float area = std::abs(tri::OutlineUtil<float>::Outline2Area(normalizedOutline));
    float qualityIndex = (perimeter * perimeter) / (4.0 * M_PI * area);
    return qualityIndex;
}

