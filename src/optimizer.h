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

#include "uv.h"

bool ParameterizeChartFromInitialTexCoord(Mesh &m, GraphManager::ChartHandle ch)
{
    auto CCIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(m, CCIDh));

    auto ICCh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "InitialConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(m, ICCh));

    PMesh pm;
    // Count the required vertices, split those at a seam in the initial parameterization
    std::unordered_map<Mesh::VertexPointer,std::unordered_map<RegionID,PMesh::VertexPointer>> mv_to_pmv{ch->FN() * 3};

    std::size_t vn = 0;
    for (auto fptr : ch->fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (mv_to_pmv.count(fptr->V(i)) == 0 || mv_to_pmv[fptr->V(i)].count(ICCh[fptr]) == 0) {
                vn++;
                mv_to_pmv[fptr->V(i)].insert(std::make_pair(ICCh[fptr], nullptr));
            }
        }
    }

    auto pmvi = tri::Allocator<PMesh>::AddVertices(pm, vn);
    auto pmfi = tri::Allocator<PMesh>::AddFaces(pm, ch->FN());

    for (auto fptr : ch->fpVec) {
        PMesh::FacePointer pmf = &*pmfi++;
        for (int i = 0; i < 3; ++i) {
            Mesh::VertexPointer mv = fptr->V(i);
            PMesh::VertexPointer& pmv = mv_to_pmv[mv][ICCh[fptr]];
            if (pmv == nullptr) pmv = &*pmvi++;
            pmf->V(i) = pmv;

            // TODO maybe add a parameter to choose the parameterization source (mesh geometry or old parameterization)
            pmv->P() = PMesh::VertexType::CoordType{fptr->cWT(i).U(), fptr->cWT(i).V(), 0};
            //pmv->P() = fptr->V(i)->P();
        }
    }

    for (auto fptr : ch->fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (face::IsBorder(*fptr, i) || CCIDh[fptr] != CCIDh[fptr->FFp(i)]) { // Border in texture space
                for (auto& idVertPair : mv_to_pmv[fptr->V0(i)]) idVertPair.second->SetB();
                for (auto& idVertPair : mv_to_pmv[fptr->V1(i)]) idVertPair.second->SetB();
            }
        }
    }

    DCPSolver<PMesh> solver(pm);
    for (auto& entry : mv_to_pmv) {
        if (entry.second.size() > 1) {
            std::vector<PMesh::VertexPointer> vertGroup;
            for (auto& idVertPair : entry.second) {
                vertGroup.push_back(idVertPair.second);
            }
            solver.DeclareUnique(vertGroup);
        }
    }

    // TEST CODE
    bool solved = solver.Solve(DCPSolver<PMesh>::__trust_vertex_border_flags);
    //bool solved = solver.Solve();

    if (solved) { // Copy texture coords back
        for (auto& fptr : ch->fpVec) {
            for (int i = 0; i < 3; ++i) {
                auto vi = fptr->V(i);
                assert(mv_to_pmv[vi].size() > 0);
                auto pmv = mv_to_pmv[vi].begin()->second;
                fptr->WT(i).P() = pmv->T().P();
            }

        }
    }

#ifndef NDEBUG

    if (!solved) std::cout << "Not solved" << std::endl;

    int n1 = 0;
    int n2 = 0;
    int n3 = 0;

    for (auto &f : pm.face) {
        Triangle3<float> tri {
            Point3f { f.WT(0).P().X(), f.WT(0).P().Y(), 0.0f },
            Point3f { f.WT(1).P().X(), f.WT(1).P().Y(), 0.0f },
            Point3f { f.WT(2).P().X(), f.WT(2).P().Y(), 0.0f }
        };
        float area = vcg::DoubleArea(tri);
        if (area < 0.0f) n1++;
        else if (area > 0.0f) n2++;
        else n3++;
    }

    std::cout << "Negative area " << n1 << std::endl;
    std::cout << "Positive area " << n2 << std::endl;
    std::cout << "Zero area     " << n3 << std::endl;

    //tri::io::ExporterOBJ<PMesh>::Save(pm, "test_mesh.obj", tri::io::Mask::IOM_WEDGTEXCOORD);

#endif

    return solved;
}

// Computes the texture outlines of a given chart
void ChartOutlinesUV(Mesh& m, GraphManager::ChartHandle chart, std::vector<std::vector<Point2f>> &outline2Vec)
{
    struct FaceAdjacency {
        bool changed[3] = {false, false, false};
        Mesh::FacePointer FFp[3];
        char FFi[3];
    };

    auto CCIDh = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(m, CCIDh));

    std::unordered_map<Mesh::FacePointer,FaceAdjacency> savedAdj;
    for (auto fptr : chart->fpVec) {
        fptr->ClearV();
        for (int i = 0; i < 3; ++i) {
            if (CCIDh[fptr] != CCIDh[fptr->FFp(i)]) {
                savedAdj[fptr].changed[i] = true;
                savedAdj[fptr].FFp[i] = fptr->FFp(i);
                savedAdj[fptr].FFi[i] = fptr->FFi(i);
                fptr->FFp(i) = fptr;
                fptr->FFi(i) = i;
            }
        }
    }

    outline2Vec.clear();
    std::vector<Point2f> outline;

    for (auto fptr : chart->fpVec) {
        for (int i = 0; i < 3; ++i) {
            if (!fptr->IsV() && face::IsBorder(*fptr, i)) {
                face::Pos<Mesh::FaceType> p(fptr, i);
                face::Pos<Mesh::FaceType> startPos = p;
                assert(p.IsBorder());
                do {
                    assert(p.IsManifold());
                    p.F()->SetV();
                    //outline.push_back(Point2<ScalarType>(p.V()->P()));
                    outline.push_back(p.F()->WT(p.VInd()).P());
                    p.NextB();
                }
                while (p != startPos);
                outline2Vec.push_back(outline);
                outline.clear();
            }
        }
    }

    for (auto& entry : savedAdj) {
        for (int i = 0; i < 3; ++i) {
            if (entry.second.changed[i]) {
                entry.first->FFp(i) = entry.second.FFp[i];
                entry.first->FFi(i) = entry.second.FFi[i];
            }
        }
    }
}

vcg::Box2f UVBox(GraphManager::ChartHandle chart) {
    vcg::Box2f bbox;
    for (auto fptr : chart->fpVec) {
        for (int i = 0; i < 3; ++i) {
            bbox.Add(fptr->cWT(i).P());
        }
    }
    return bbox;
}

// returns the number of charts that could not be parameterized
/// TODO update distortion info if needed (this should also be done through the graph manager)
static int ParameterizeGraph(std::shared_ptr<MeshGraph> graph)
{
    Timer timer;
    Mesh& m = graph->mesh;

    std::cout << "Parameterizing " << graph->charts.size() << " regions..." << std::endl;

    std::vector<std::vector<Point2f>> texOutlines;
    std::unordered_map<RegionID,std::size_t> outlineMap; // map each region to the index of its outline in texOutlines
    int iter = 0;
    int numFailed = 0;
    for (auto entry : graph->charts) {
        std::shared_ptr<FaceGroup> chart = entry.second;

        chart->NotifyParameterizationChange(); // always changes due to repacking, even if the original coordinates were preserved

        std::cout << "Chart " << chart->id << " - FN=" << chart->FN() << ", FI=" << tri::Index(m, chart->Fp()) << std::endl;

        // If the old chart has no valid parameterization, skip it
        float oldUvArea = chart->AreaUV();
        if (oldUvArea == 0) continue;

        if (chart->numMerges == 0) {
            // The chart was not merged, therefore there is nothing to do -- the original parameterization is preserved

            // !!! this does not happen to work as there are issues with selecting the outline2 from such regions
            // Example: closed surfaces parameterized without cuts (like two sides projected to a plane)

            // TODO a possible workaround would be to use the bounding box as the outline
        }
        else {
            // TODO function parameter to choose the source of the mesh geometry (3D or old UV space)
            bool parameterized = ParameterizeChartFromInitialTexCoord(m, chart);

            if (parameterized) {
                // Normalize area: the region gets scaled by sqrt(oldUVArea)/sqrt(newUVArea) to keep the original proportions
                // between the regions of the atlas in an attempt to not lose too much texture data when the new texture is rendered
                float newUvArea = chart->AreaUV();
                float scale = std::sqrt(oldUvArea / newUvArea);
                assert(scale > 0);
                vcg::Box2f uvBox = UVBox(chart);
                for (auto fptr : chart->fpVec) {
                    for (int i = 0; i < 3; ++i) {
                        fptr->WT(i).P() = (fptr->WT(i).P() - uvBox.min) * scale;
                    }
                }

                chart->numMerges = 0;
            }
            else {
                numFailed++;
                std::cout << "Parameterization failed" << std::endl;
            }
        }
        // Save the outline of the parameterization for this portion of the mesh

        std::vector<std::vector<Point2f>> uvOutlines;
        ChartOutlinesUV(m, chart, uvOutlines);
        int i = tri::OutlineUtil<float>::LargestOutline2(uvOutlines);
        if (tri::OutlineUtil<float>::Outline2Area(uvOutlines[i]) < 0)
            tri::OutlineUtil<float>::ReverseOutline2(uvOutlines[i]);

        outlineMap[chart->id] = texOutlines.size();
        texOutlines.push_back(uvOutlines[i]);

        std::cout << "Iteration " << iter++ << " took " << timer.TimeSinceLastCheck() << " seconds" << std::endl;
    }

    std::cout << "Parameterization took " << timer.TimeElapsed() << " seconds" << std::endl;

    std::cout << "Packing the atlas..." << std::endl;

    // pack the atlas TODO function parameter to choose the packing strategy
    RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters packingParam;
    packingParam.costFunction  = RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters::LowestHorizon;
    packingParam.doubleHorizon = true;
    packingParam.cellSize = 4;
    packingParam.rotationNum = 16; //number of rasterizations in 90°

    Point2i gridSize(1024, 1024);
    std::vector<Similarity2f> transforms;

    //RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Pack(texOutlines, gridSize, transforms, packingParam);

    Point2f cover;
    //PolyPacker<float>::PackAsObjectOrientedRect(texOutlines, gridSize, transforms, cover);
    PolyPacker<float>::PackAsAxisAlignedRect(texOutlines, gridSize, transforms, cover);

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

    return numFailed;
}

/// TODO refactoring
static void ReduceTextureFragmentation_NoPacking(GraphManager &gm, std::size_t minRegionSize);

static void ReduceTextureFragmentation(Mesh &m, std::shared_ptr<MeshGraph> graph, std::size_t minRegionSize)
{
    if (minRegionSize == 0) return;

    assert(minRegionSize < (std::size_t) m.FN());

    tri::UpdateTopology<Mesh>::FaceFace(m);

    GraphManager gm{graph};

    std::cout << "Initialized graph manager" << std::endl;

    ReduceTextureFragmentation_NoPacking(gm, minRegionSize);

    int c = ParameterizeGraph(gm.Graph());
    if (c > 0) std::cout << "WARNING: " << c << " regions were not parameterized correctly" << std::endl;
}

static void ReduceTextureFragmentation_NoPacking(GraphManager &gm, std::size_t minRegionSize)
{
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
}

static void ReduceTextureFragmentation2(Mesh &m, std::shared_ptr<MeshGraph> graph, std::size_t minRegionSize)
{
    if (minRegionSize == 0) return;

    assert(minRegionSize < (std::size_t) m.FN());

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
        Timer t_;
        PMesh pm;
        for (auto fptr : chart->fpVec) {
            tri::Allocator<PMesh>::AddFace(pm, fptr->P(0), fptr->P(1), fptr->P(2));
        }
        tri::Clean<PMesh>::RemoveDuplicateVertex(pm);
        tri::Allocator<PMesh>::CompactEveryVector(pm);
        std::cout << "Per face mesh alloc took " << t_.TimeElapsed() << std::endl;
        //tri::UpdateTopology<PMesh>::FaceFace(pm);
        tri::UpdateSelection<PMesh>::Clear(pm);

        //tri::io::ExporterOBJ<PMesh>::Save(pm, "submesh_iteration.obj", tri::io::Mask::IOM_WEDGTEXCOORD);

        tri::PoissonSolver<PMesh> solver(pm);
        //DCPSolver<PMesh> solver(pm);
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

