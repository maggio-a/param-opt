#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/poisson_solver.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <vcg/space/rasterized_outline2_packer.h>
#include <vcg/space/outline2_packer.h>
#include <wrap/qt/outline2_rasterizer.h>

#include <vector>

#include "mesh_graph.h"
#include "timer.h"
#include "dcpsolver.h"

static void ReduceTextureFragmentation(Mesh &m, std::shared_ptr<MeshGraph> graph, std::size_t minRegionSize)
{
    if (minRegionSize == 0) return;

    assert(minRegionSize < m.FN());

    tri::UpdateTopology<Mesh>::FaceFace(m);

    GraphManager gm{graph};

    std::cout << "Initialized graph manager" << std::endl;

    Timer timer;
    int mergeCount;
    int numIter = 0;

    do {
        mergeCount = gm.CloseMacroRegions(minRegionSize);

        while (gm.HasNextEdge()) {
            auto we = gm.PeekNextEdge();
            if (we.first.a->FN() > minRegionSize && we.first.b->FN() > minRegionSize)
                break;
            else {
                gm.RemoveNextEdge();
                gm.Collapse(we.first);
                mergeCount++;
                if (mergeCount%50 == 0) {
                    std::cout << "Merged " << mergeCount << " regions..." << std::endl;
                }
            }
        }

        std::cout << "Iteration "  << numIter << " took " << timer.TimeSinceLastCheck() << " seconds ("
                  << mergeCount << " merges took place)" << std::endl;

        numIter++;

    } while (mergeCount > 0);

    std::cout << "Stopping after " << numIter << " passes and " << timer.TimeElapsed() << " seconds" << std::endl;

    timer.Reset();

    std::cout << "Parameterizing " << graph->charts.size() << " regions..." << std::endl;
    int iter = 0;

    std::vector<std::vector<Point2f>> texOutlines;
    std::unordered_map<std::size_t, std::size_t> outlineMap;
    // reparameterize each region with the poisson solver
    for (auto entry : graph->charts) {
        std::shared_ptr<FaceGroup> chart = entry.second;

        std::cout << "Chart " << chart->id << " - FN=" << chart->FN() << ", FI=" << tri::Index(m, chart->Fp()) << std::endl;

        // If the old chart has no valid parameterization, skip it
        float oldUvArea = chart->AreaUV();
        if (oldUvArea == 0) continue;

        // Build the mesh to parameterize
        PMesh pm;
        for (auto fptr : chart->fpVec) {
            tri::Allocator<PMesh>::AddFace(pm, fptr->P(0), fptr->P(1), fptr->P(2));
        }
        tri::Clean<PMesh>::RemoveDuplicateVertex(pm);
        tri::Allocator<PMesh>::CompactEveryVector(pm);
        //tri::UpdateTopology<PMesh>::FaceFace(pm);
        tri::UpdateSelection<PMesh>::Clear(pm);

        //tri::io::ExporterOBJ<PMesh>::Save(pm, "submesh_iteration.obj", tri::io::Mask::IOM_WEDGTEXCOORD);

        tri::PoissonSolver<PMesh> solver(pm);
        //DCPSolver<PMesh> solver(pm);//  assert scale < 0
        tri::UpdateBounding<PMesh>::Box(pm);
        if (!solver.IsFeasible()) {

            // Note: since we merge two regions only if the result is parameterizable
            // with the PoissonSolver, unfeasibility can only occur for regions that have NOT been merged with anything
            // else, thus having the original parameterization untouched. Just assert the condition and copy the original
            // parameterization over to the PMesh

            assert(chart->numMerges == 0);

            for (std::size_t i = 0; i < pm.face.size(); ++i) {
                auto &pf = pm.face[i];
                auto &f = *(chart->fpVec[i]);
                for (int k = 0; k < f.VN(); ++k) {
                    pf.WT(k) = f.WT(k);
                }
            }

            // !!! this does not happen to work as there are issues with selecting the outline2 from such regions
            // Example: closed surfaces parameterized without cuts (like two sides projected to a plane)

        } else {

            // In this case generate a new parameterization

            // Solve the system and update texture coords
            solver.Init();
            solver.FixDefaultVertices();
            solver.SolvePoisson();
            //tri::UpdateTexture<PMesh>::WedgeTexFromVertexTex(pm);

            // Normalize area: the region gets scaled by sqrt(oldUVArea)/sqrt(newUVArea) to keep the original proportions
            // between the regions of the atlas in an attempt to not lose too much texture data when the new texture is rendered
            float newUvArea = 0.0f;
            for (auto &pf : pm.face) {
                newUvArea += std::abs(tri::Distortion<PMesh,true>::AreaUV(&pf));
            }

            float scale = std::sqrt(oldUvArea / newUvArea);

            /*if (scale <= 0) {
                char name[50] = {0};
                std::sprintf(name, "chart_%d.obj", iter);
                tri::io::ExporterOBJ<PMesh>::Save(pm, name, tri::io::Mask::IOM_WEDGTEXCOORD);
            }*/

            assert(scale > 0);

            vcg::Box2f uvBox = tri::UV_Utils<PMesh>::PerWedgeUVBox(pm);
            for (std::size_t i = 0; i < pm.face.size(); ++i) {
                auto &pf = pm.face[i];
                auto &f = *(chart->fpVec[i]);
                for (int k = 0; k < f.VN(); ++k) {
                    pf.WT(k).P() = (pf.WT(k).P() - uvBox.min) * scale;
                    f.WT(k) = pf.WT(k);
                }
            }
        }

        // Save the outline of the parameterization for this portion of the mesh

        std::vector<std::vector<Point2f>> regionOutlines;
        ConvertTextureBoundaryToOutline2Vec<float, PMesh>(pm, regionOutlines);
        int i = tri::OutlineUtil<float>::LargestOutline2(regionOutlines);
        if (tri::OutlineUtil<float>::Outline2Area(regionOutlines[i]) < 0)
            tri::OutlineUtil<float>::ReverseOutline2(regionOutlines[i]);

        outlineMap[entry.first] = texOutlines.size();
        texOutlines.push_back(regionOutlines[i]);

        std::cout << "Iteration " << iter++ << " took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    }

    std::cout << "Parameterization took " << timer.TimeElapsed() << " seconds" << std::endl;

    std::cout << "Packing the atlas..." << std::endl;

    // pack the atlas
    RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters packingParam;
    packingParam.costFunction  = RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters::LowestHorizon;
    packingParam.doubleHorizon = true;
    packingParam.cellSize = 4;
    packingParam.rotationNum = 16; //number of rasterizations in 90°

    Point2i gridSize(1024, 1024);
    std::vector<Similarity2f> transforms;

    RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Pack(texOutlines, gridSize, transforms, packingParam);

    Point2f cover;
    //PolyPacker<float>::PackAsObjectOrientedRect(texOutlines, gridSize, transforms, cover);
    //PolyPacker<float>::PackAsAxisAlignedRect(texOutlines, gridSize, transforms, cover);

    std::cout << "Packing took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;

    //assert(transforms.size() == pdata.charts.size());

    for (auto p : outlineMap) {
        for (auto fptr : graph->charts[p.first]->fpVec) {
            for (int j = 0; j < fptr->VN(); ++j) {
                Point2f transformedTexCoordPos = transforms[p.second] * (fptr->WT(j).P());
                fptr->WT(j).P() = transformedTexCoordPos / 1024.0f;
            }
        }
    }

}

#endif // OPTIMIZER_H

