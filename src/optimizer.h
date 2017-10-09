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

using namespace vcg;

using ChartMap = std::unordered_map<size_t, shared_ptr<FaceGroup>>;

static void MergeRegions(Mesh &m, std::shared_ptr<FaceGroup> pc1, std::shared_ptr<FaceGroup> pc2, ChartMap& regions)
{
    auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(m, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(m, CHIDh));

    // update id and add face to other region
    for (auto fp : pc2->fpVec) {
        CHIDh[fp] = pc1->id;
        pc1->AddFace(fp, CHIDh);
    }

    // update adjacencies
    pc1->adj.erase(pc2);
    for (auto regionptr : pc2->adj) {
        if (regionptr->id != pc1->id) {
            regionptr->adj.erase(pc2);
            regionptr->adj.insert(pc1);
            pc1->adj.insert(regionptr);
        }
    }

    regions.erase(pc2->id);

    pc1->numMerges++;
}

struct MergeStep {
    bool feasible;
    std::shared_ptr<FaceGroup> baseRegion;
    float quality;
};

/// Builds a PMesh with face face topology and bounding box initialized, so that it can be
/// passed to the poisson solver
template <typename MeshType>
static void BuildPMeshFromFacePointers(PMesh &pm, const std::vector<std::vector<typename MeshType::FacePointer>* >& vFpVecp)
{
    pm.Clear();

    auto f = [&pm](typename MeshType::FacePointer fptr) {
        tri::Allocator<PMesh>::AddFace(pm, fptr->P(0), fptr->P(1), fptr->P(2));
    };

    for (auto fpVecp : vFpVecp) std::for_each(fpVecp->begin(), fpVecp->end(), f);

    tri::Clean<PMesh>::RemoveDuplicateVertex(pm);
    tri::Allocator<PMesh>::CompactEveryVector(pm);

    tri::UpdateTopology<PMesh>::FaceFace(pm);
    tri::UpdateBounding<PMesh>::Box(pm);
}

/// Returns true if the mesh can be parameterized by the poisson solver
template <typename MeshType>
static bool Parameterizable(MeshType &m)
{
    if (tri::Clean<MeshType>::CountNonManifoldEdgeFF(m) > 0) {
        return false;
    }

    if (tri::Clean<MeshType>::CountNonManifoldVertexFF(m) > 0) {
        return false;
    }

    if (tri::Clean<MeshType>::IsWaterTight(m)) {
        return false;
    }

    if (tri::Clean<MeshType>::MeshGenus(m) > 0) {
        return false;
    }

    return true;
}

// Evaluates the possible merge steps, returns true if at least one is feasible
static bool EvaluateMergeSteps(Mesh &m, std::shared_ptr<FaceGroup> chart, std::vector<MergeStep> &steps)
{
    steps.clear();
    steps.reserve(chart->NumAdj());

    bool feasible = false;
    assert(tri::HasPerFaceAttribute(m, "ConnectedComponentID"));

    // First evaluate shared border
    std::unordered_map<RegionID, float> borderMap;
    auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<std::size_t>(m, "ConnectedComponentID");
    for (auto fptr : chart->fpVec) {
        for (int i = 0; i < fptr->VN(); ++i) {
            auto ffpi = fptr->FFp(i);
            auto adjId = CHIDh[ffpi];
            if (adjId != CHIDh[fptr]) {
                borderMap[adjId] += DistortionWedge::EdgeLenght3D(fptr, i);
            }
        }
    }

    // Then gather data for each possible merge step
    for (auto &adjacentRegion : chart->adj) {
        assert(borderMap.count(adjacentRegion->id));

        PMesh probe;
        BuildPMeshFromFacePointers<Mesh>(probe, {&(chart->fpVec), &(adjacentRegion->fpVec)});

        bool feasibleMove = Parameterizable<PMesh>(probe);

        if (feasibleMove) {
            feasible = true;
            tri::PoissonSolver<PMesh> ps(probe);
            ps.Init();
            ps.FixDefaultVertices();
            bool solved = ps.SolvePoisson();
            assert(solved);

            float areaScaleFactor, edgeScaleFactor;
            float areaDistortion = 0;
            tri::Distortion<PMesh,true>::MeshScalingFactor(probe, areaScaleFactor, edgeScaleFactor);
            for (auto& face : probe.face) {
                areaDistortion += tri::Distortion<PMesh,true>::AreaDistortion(&face, areaScaleFactor);
            }
            steps.push_back(MergeStep{true, adjacentRegion, 1.0f / (areaDistortion / probe.FN())});
        } else {
            steps.push_back(MergeStep{false, nullptr, 0});
        }

    }

    return feasible;
}

/// Heuristic procedure that simulates a 'close' operation over large texture regions
/// in order to remove multiple islands in one go. This can be useful if there are stalemates
/// where the presence of multiple islands makes the region outline non manifold. In such
/// a case merging the islands together rather than one at a time can fix the issue
static int CloseMacroRegions(Mesh& m, ChartMap& regions, std::size_t minRegionSize)
{
    int mergeCount = 0;

    // Operate on two passes, first build an index of the merges and an inverted index of the already merged small regions
    // then perform the actual merges

    std::unordered_map<RegionID, std::vector<RegionID>> mergeLists;
    std::unordered_map<RegionID, RegionID> invertedIndex;

    for (const auto& entry : regions) {
        auto chart = entry.second;
        if (invertedIndex.count(chart->id) == 1) continue; // skip if this region is already going to be merged to something else
        for (auto& adjRegion : chart->adj) if (adjRegion->NumAdj() == 1 && adjRegion->FN() < minRegionSize) {
            assert(invertedIndex.count(adjRegion->id) == 0);
            mergeLists[chart->id].push_back(adjRegion->id);
            invertedIndex[adjRegion->id] = chart->id;
        }
    }

    for (auto& entry : mergeLists) {
        PMesh probe;
        std::vector<std::vector<Mesh::FacePointer>* > fpVecp;
        fpVecp.reserve(entry.second.size()+1);
        fpVecp.push_back(&(regions[entry.first]->fpVec));
        for (auto id : entry.second) {
            fpVecp.push_back(&(regions[id]->fpVec));
        }
        BuildPMeshFromFacePointers<Mesh>(probe, fpVecp);
        if (Parameterizable(probe)) {
            std::cout << "Merging " << entry.second.size() << " islands from macro region " << entry.first << std::endl;
            for (auto id : entry.second) {
                MergeRegions(m, regions[entry.first], regions[id], regions);
                mergeCount++;
            }
        } else {
            //tri::io::ExporterOBJ<PMesh>::Save(probe, "fail.obj", tri::io::Mask::IOM_WEDGTEXCOORD);
        }
    }

    return mergeCount;
}

static void ReduceTextureFragmentation(Mesh &m, MeshGraph &pdata, std::size_t minRegionSize)
{
    using ChartMap = decltype(pdata.charts);
    using CVal = ChartMap::value_type;
    using BVal = std::unordered_map<std::size_t, float>::value_type;

    if (minRegionSize == 0) return;

    assert(minRegionSize < m.FN());

    tri::UpdateTopology<Mesh>::FaceFace(m);

    Timer timer;
    int mergeCount;
    int numIter = 0;

    do {
        mergeCount = CloseMacroRegions(m, pdata.charts, minRegionSize);

        // Filter regions that are not parameterized or that have no adjacencies
        ChartMap regions;
        ChartMap savedRegions; // TODO rename to discardedRegions

        std::partition_copy(pdata.charts.begin(), pdata.charts.end(),
                            std::inserter(regions, regions.end()), std::inserter(savedRegions, savedRegions.end()),
                            [minRegionSize](const CVal &val) { return val.second->FN() < minRegionSize && val.second->NumAdj() > 0; });

        while (regions.size() > 0) {

            auto chartEntry = regions.begin();
            assert(chartEntry != regions.end());

            bool merged = false;

            auto chart = chartEntry->second;

            if (chart->FN() < minRegionSize) {
                // The region is a candidate to be merged
                std::vector<MergeStep> steps;
                bool validSteps = EvaluateMergeSteps(m, chart, steps);
                if (validSteps) {
                    // Sort possible merge steps in descending order of quality
                    std::sort(steps.begin(), steps.end(), [](const MergeStep& m1, const MergeStep& m2) { return m1.quality > m2.quality; });
                    auto &step = steps[0];
                    assert(step.feasible);
                    MergeRegions(m, step.baseRegion, chart, regions);

                    mergeCount++;
                    merged = true;

                    if (mergeCount%50 == 0) {
                        std::cout << "Merged " << mergeCount << " regions..." << std::endl;
                    }
                }
            }

            if (!merged) {
                // In this case the current region cannot be merged to any neighboring one, so remove it
                savedRegions[chartEntry->first] = regions[chartEntry->first];
                regions.erase(chartEntry->first);
            }
        } // while (regions.size() > 0)

        std::copy(savedRegions.begin(), savedRegions.end(), std::inserter(regions, regions.end()));

        std::cout << "Performing final closing pass..." << std::endl;
        mergeCount += CloseMacroRegions(m, regions, minRegionSize);

        numIter++;

        std::cout << "Iteration "  << numIter << " took " << timer.TimeSinceLastCheck() << " seconds ("
                  << mergeCount << " merges took place)" << std::endl;

        pdata.charts = regions;

    } while (mergeCount > 0);

    std::cout << "Stopping after " << numIter << " passes and " << timer.TimeElapsed() << " seconds" << std::endl;

    timer.Reset();

    std::cout << "Parameterizing " << pdata.charts.size() << " regions..." << std::endl;
    int iter = 0;

    std::vector<std::vector<Point2f>> texOutlines;
    std::unordered_map<std::size_t, std::size_t> outlineMap;
    // reparameterize each region with the poisson solver
    for (auto entry : pdata.charts) {
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

        } else {

            // In this case generate a new parameterization

            // Solve the system and update texture coords
            solver.Init();
            solver.FixDefaultVertices();
            solver.SolvePoisson();
            tri::UpdateTexture<PMesh>::WedgeTexFromVertexTex(pm);

            // Normalize area: the region gets scaled by sqrt(oldUVArea)/sqrt(newUVArea) to keep the original proportions
            // between the regions of the atlas in an attempt to not lose too much texture data when the new texture is rendered
            float newUvArea = 0.0f;
            for (auto &pf : pm.face) {
                newUvArea += std::abs(tri::Distortion<PMesh,true>::AreaUV(&pf));
            }

            float scale = std::sqrt(oldUvArea / newUvArea);

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
        for (auto fptr : pdata.charts[p.first]->fpVec) {
            for (int j = 0; j < fptr->VN(); ++j) {
                Point2f transformedTexCoordPos = transforms[p.second] * (fptr->WT(j).P());
                fptr->WT(j).P() = transformedTexCoordPos / 1024.0f;
            }
        }
    }

}

#endif // OPTIMIZER_H

