#include "mesh.h"
#include "mesh_graph.h"
#include "uv.h"
#include "texture_optimization.h"
#include "texture_rendering.h"
#include "timer.h"
#include "gl_utils.h"
#include "texture.h"
#include "packing_utils.h"
#include "logging.h"

#include <wrap/io_trimesh/io_mask.h>

#include <string>
#include <vector>
#include <iostream>

#include <QCoreApplication>
#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>


using namespace vcg;

int main(int argc, char *argv[])
{
    // Make sure the executable directory is added to Qt's library path
    {
        QFileInfo executableInfo(argv[0]);
        QCoreApplication::addLibraryPath(executableInfo.dir().absolutePath());
    }

    Mesh m;
    TextureObjectHandle textureObject;
    int loadMask;

    ensure_condition(argc > 1);

    if (LoadMesh(m, argv[1], textureObject, loadMask) == false) {
        std::cout << "Failed to open mesh" << std::endl;
        std::exit(-1);
    }

    ensure_condition(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    tri::UpdateTopology<Mesh>::FaceFace(m);

    auto graph = ComputeParameterizationGraph(m, textureObject);


    ScaleTextureCoordinatesToImage(m, textureObject);
    StoreWedgeTexCoordAsAttribute(m);
    ScaleTextureCoordinatesToImage(m, textureObject);


    std::vector<RegionID> nullRegions;
    for (auto entry : graph->charts) {
        ChartHandle chart = entry.second;
        if (chart->AreaUV() == 0) {
            nullRegions.push_back(chart->id);
        }
    }
    for (auto id : nullRegions) {
        graph->charts.erase(id);
    }


    GLInit();

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
        ScaleTextureCoordinatesToImage(m, textureObject);

        PackingOptions opt = { cost[idx], false, false, false, false };
        stats.push_back(Pack(graph, opt));
        std::vector<RasterizedParameterizationStats> after = GetRasterizationStats(m, textureObject);
        std::cout << "[LOG] Raster stats after parameterizing" << std::endl;
        LogParameterizationStats(graph, after);
        TextureObjectHandle newTexture = RenderTexture(m, textureObject, false, InterpolationMode::Linear, nullptr);
        std::string savename = "out_" + costName[idx] + std::string("_hi") + m.name;
        if (SaveMesh(m, savename.c_str(), newTexture, true) == false) {
            std::cout << "Model not saved correctly" << std::endl;
        }
    }

    {
        ScaleTextureCoordinatesToImage(m, textureObject);

        PackingOptions opt = { cost[idx], false, false, true, false };
        stats.push_back(Pack(graph, opt));
        std::vector<RasterizedParameterizationStats> after = GetRasterizationStats(m, textureObject);
        std::cout << "[LOG] Raster stats after parameterizing" << std::endl;
        LogParameterizationStats(graph, after);
        TextureObjectHandle newTexture = RenderTexture(m, textureObject, false, InterpolationMode::Linear, nullptr);
        std::string savename = "out_" + costName[idx] + std::string("_hi_inner") + m.name;
        if (SaveMesh(m, savename.c_str(), newTexture, true) == false) {
            std::cout << "Model not saved correctly" << std::endl;
        }
    }

    {
        ScaleTextureCoordinatesToImage(m, textureObject);

        PackingOptions opt = { cost[idx], true, true, true, false };
        stats.push_back(Pack(graph, opt));
        std::vector<RasterizedParameterizationStats> after = GetRasterizationStats(m, textureObject);
        std::cout << "[LOG] Raster stats after parameterizing" << std::endl;
        LogParameterizationStats(graph, after);
        TextureObjectHandle newTexture = RenderTexture(m, textureObject, false, InterpolationMode::Linear, nullptr);
        std::string savename = "out_" + costName[idx] + std::string("_low_perm_inner") + m.name;
        if (SaveMesh(m, savename.c_str(), newTexture, true) == false) {
            std::cout << "Model not saved correctly" << std::endl;
        }
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
