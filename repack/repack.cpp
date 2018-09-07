#include "mesh.h"
#include "mesh_graph.h"
#include "uv.h"
#include "texture_optimization.h"
#include "texture_rendering.h"
#include "parameterization_checker.h"
#include "timer.h"
#include "gl_utils.h"
#include "texture.h"
#include "packing_utils.h"

#include <wrap/io_trimesh/io_mask.h>

#include <string>
#include <vector>
#include <iostream>

#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>


using namespace vcg;


void LogParameterizationStats(std::shared_ptr<MeshGraph> graph, RasterizedParameterizationStats stats, const std::string& header)
{
    int boundaries;
    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(graph->mesh);
    boundaries = tri::Clean<Mesh>::CountHoles(graph->mesh);
    tri::UpdateTopology<Mesh>::FaceFace(graph->mesh);

    std::cout << header << std::endl;
    std::cout << "[LOG] Texture Size , Charts , Boundaries , Occupancy , BorderUV(px), AreaUV(px), BilinearAreaUV(px), Overwritten fragments, Lost Fragments" << std::endl;
    std::cout << "[LOG] "
              << stats.rw << "x" << stats.rh << " , "
              << graph->charts.size() << " , "
              << boundaries << " , "
              << (stats.totalFragments - stats.lostFragments) / double(stats.rw*stats.rh) << " , "
              << stats.boundaryFragments << " , "
              << stats.totalFragments  << " , "
              << stats.totalFragments_bilinear  << " , "
              << stats.overwrittenFragments  << " , "
              << stats.lostFragments << "" << std::endl;

}

void CorrectAspectRatio(Mesh& m, TextureObjectHandle textureObject)
{
    std::vector<double> ratios;
    for (std::size_t i = 0; i < textureObject->ArraySize(); ++i) {
        double uvRatio = textureObject->TextureWidth(i) / (double) textureObject->TextureHeight(i);
        ratios.push_back(uvRatio);
    }

    for (auto &f : m.face) {
        int textureIndex = f.WT(0).N();
        double uvRatio = ratios[textureIndex];
        for (int i = 0; i < 3; ++i) {
            if (uvRatio > 1)
                f.WT(i).P().X() *= uvRatio;
            else if (uvRatio < 1)
                f.WT(i).P().Y() /= uvRatio;
        }
    }
}

int main(int argc, char *argv[])
{
    Mesh m;
    TextureObjectHandle textureObject;
    int loadMask;

    assert(argc > 1);

    if (LoadMesh(m, argv[1], textureObject, loadMask) == false) {
        std::cout << "Failed to open mesh" << std::endl;
        std::exit(-1);
    }

    assert(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    tri::UpdateTopology<Mesh>::FaceFace(m);

    auto graph = ComputeParameterizationGraph(m, textureObject);

    StoreWedgeTexCoordAsAttribute(m);

    PrintParameterizationInfo(graph);

    GLInit();

    CorrectAspectRatio(m, textureObject);

    /*
    std::unordered_map<RegionID,float> qualityMap;
    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(m);
    for (auto& entry : graph->charts) {
        Timer t;
        std::vector<std::vector<Point2f>> uvOutlines;
        ChartOutlinesUV(graph->mesh, *(entry.second), uvOutlines);
        int i = tri::OutlineUtil<float>::LargestOutline2(uvOutlines);
        if (tri::OutlineUtil<float>::Outline2Area(uvOutlines[i]) < 0)
            tri::OutlineUtil<float>::ReverseOutline2(uvOutlines[i]);
        //qualityMap[entry.first] = EvaluatePackingQuality_Horizon(uvOutlines[i]);
        qualityMap[entry.first] = EvaluatePackingQuality_VoidArea(uvOutlines[i]);
        //qualityMap[entry.first] = EvaluatePackingQuality_Perimeter(uvOutlines[i]);
        std::cout << entry.first << " " << qualityMap[entry.first] << std::endl;
        std::cout << "Evaluation took " << t.TimeElapsed() << " seconds." << std::endl;
    }

    GeneratePackingQualityTexture(graph, textureObject, qualityMap);
    GLTerminate();
    return 0;
    */

    auto minWastedSpace = RasterizationBasedPacker::Parameters::CostFuncEnum::MinWastedSpace;
    auto lowestHorizon = RasterizationBasedPacker::Parameters::CostFuncEnum::LowestHorizon;

    RasterizationBasedPacker::Parameters::CostFuncEnum cost[] = { minWastedSpace, lowestHorizon };
    std::string costName[] = { "minWastedSpace", "lowestHorizon" };

    int idx = 0; // min wasted space
    //int idx = 1; // lowest horizon

    if (argc > 2)
        idx = atoi(argv[2]);

    std::vector<RasterizationBasedPacker::PackingStats> stats;

    {
        PackingOptions opt = { cost[idx], false, false, false, false };
        stats.push_back(Pack(graph, opt));
        RasterizedParameterizationStats after = GetRasterizationStats(m, textureObject->imgVec[0]->width(), textureObject->imgVec[0]->height());
        LogParameterizationStats(graph, after, std::string("[LOG] Raster stats after parameterizing"));
        TextureObjectHandle newTexture = RenderTexture(m, textureObject, false, InterpolationMode::Linear, nullptr);
        std::string savename = "out_" + costName[idx] + std::string("_hi") + m.name;
        if (SaveMesh(m, savename.c_str(), newTexture, true) == false) {
            std::cout << "Model not saved correctly" << std::endl;
        }

        CorrectAspectRatio(m, newTexture);
    }

    {
        PackingOptions opt = { cost[idx], false, false, true, false };
        stats.push_back(Pack(graph, opt));
        RasterizedParameterizationStats after = GetRasterizationStats(m, textureObject->imgVec[0]->width(), textureObject->imgVec[0]->height());
        LogParameterizationStats(graph, after, std::string("[LOG] Raster stats after parameterizing"));
        TextureObjectHandle newTexture = RenderTexture(m, textureObject, false, InterpolationMode::Linear, nullptr);
        std::string savename = "out_" + costName[idx] + std::string("_hi_inner") + m.name;
        if (SaveMesh(m, savename.c_str(), newTexture, true) == false) {
            std::cout << "Model not saved correctly" << std::endl;
        }

        CorrectAspectRatio(m, newTexture);
    }

    {
        PackingOptions opt = { cost[idx], true, true, true, false };
        stats.push_back(Pack(graph, opt));
        RasterizedParameterizationStats after = GetRasterizationStats(m, textureObject->imgVec[0]->width(), textureObject->imgVec[0]->height());
        LogParameterizationStats(graph, after, std::string("[LOG] Raster stats after parameterizing"));
        TextureObjectHandle newTexture = RenderTexture(m, textureObject, false, InterpolationMode::Linear, nullptr);
        std::string savename = "out_" + costName[idx] + std::string("_low_perm_inner") + m.name;
        if (SaveMesh(m, savename.c_str(), newTexture, true) == false) {
            std::cout << "Model not saved correctly" << std::endl;
        }

        CorrectAspectRatio(m, newTexture);
    }

    std::cout << "[PACKING_TABLE_HEADER], Mesh, CostFunc, HiRes2Hor, HiRes4Hor, LowResPerm4Hor, PermIdx" << std::endl;
    std::cout << "[PACKING_TABLE_ROW], " << m.name << ", " << costName[idx] << ", " << stats[0].efficiency << ", "
              << stats[1].efficiency << ", "<< stats[2].efficiency << ", "<< stats[2].permutationIndex << std::endl;

    //PrintParameterizationInfo(graph);

    //std::cout << std::endl << "  Texture occupancy difference: "
    //          << ((after.totalFragments - after.lostFragments) / double(after.rw*after.rh)) - ((before.totalFragments - before.lostFragments) / double(before.rw*before.rh)) << std::endl;

    GLTerminate();

    //auto uvBox = UVBox(m);
    //std::cout << "BC Efficiency: " << graph->AreaUV() / (double)(uvBox.DimX() * uvBox.DimY()) << std::endl;

    return 0;


}
