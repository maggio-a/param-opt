#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>
#include <vcg/complex/algorithms/parametrization/poisson_solver.h>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <string>
#include <vector>
#include <array>
#include <memory>
#include <unordered_set>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>

#include <unistd.h>

#include <QImage>

#include "mesh.h"
#include "texture_rendering.h"

using namespace vcg;

class PMeshVertex;
class PMeshFace;
struct PMeshUsedTypes : public UsedTypes<Use<PMeshVertex>::AsVertexType, Use<PMeshFace>::AsFaceType> {};

class PMeshVertex : public Vertex<PMeshUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::TexCoord2f, vertex::BitFlags> {};
class PMeshFace : public Face<PMeshUsedTypes, face::VertexRef, face::FFAdj, face::WedgeTexCoord2f, face::BitFlags> {};
class PMesh : public tri::TriMesh<std::vector<PMeshVertex>, std::vector<PMeshFace>> {};

using DistortionWedge = tri::Distortion<Mesh,true>;

inline std::string GetFileName(const std::string& path)
{
    auto pos = path.find_last_of("/");
    if (pos == std::string::npos) return path;
    else return path.substr(pos+1);
}

inline std::string GetDirName(const std::string &path)
{
    auto pos = path.find_last_of("/");
    if (pos == std::string::npos) return std::string(".");
    else return path.substr(0, pos);
}

struct ChartData {
    const std::size_t id;
    std::vector<Mesh::FacePointer> fpVec;
    std::unordered_set<std::shared_ptr<ChartData>> adj;
    float meshArea;
    float chartArea;
    float chartBorder;
    float uvMeshBorder;

    ChartData(const std::size_t id_) : id{id_}, fpVec{}, adj{}, meshArea{0}, chartArea{0}, chartBorder{0}, uvMeshBorder{0} {}

    void AddFace(const Mesh::FacePointer fptr, Mesh::PerFaceAttributeHandle<std::size_t>& CCIDh) {
        fpVec.push_back(fptr);
        meshArea += tri::Distortion<Mesh,true>::Area3D(fptr);
        chartArea += std::abs(tri::Distortion<Mesh,true>::AreaUV(fptr));
        int vn = fptr->VN();
        for (int i = 0; i < vn; ++i) {
            auto ffpi = fptr->FFp(i);
            if (id != CCIDh[ffpi])
                chartBorder += (fptr->cWT((i+1)%vn).P() - fptr->cWT(i).P()).Norm();
            else if (ffpi == fptr)
                uvMeshBorder += (fptr->cWT((i+1)%vn).P() - fptr->cWT(i).P()).Norm();
        }
    }

    float ChartArea() {
        float areaUV = 0.0f;
        for (auto fp : fpVec) areaUV += tri::Distortion<Mesh,true>::AreaUV(fp);
        return areaUV;
    }

    Mesh::FacePointer Fp() {
        assert(!fpVec.empty());
        return fpVec[0];
    }

    std::size_t FN() {
        return fpVec.size();
    }
};

struct ParameterizationData {

    std::unordered_map<std::size_t, std::shared_ptr<ChartData>> charts;
    std::vector<std::shared_ptr<QImage>> textures;

    ParameterizationData() : charts{}, textures{} {}

    std::shared_ptr<ChartData> GetChart(std::size_t i)
    {
        if (charts.find(i) == charts.end()) charts.insert(std::make_pair(i, std::make_shared<ChartData>(i)));
        return charts[i];
    }

    std::size_t Count() {
        std::size_t sz = 0;
        for (auto c : charts) {
            sz += c.second->fpVec.size();
        }
        return sz;
    }

    float Area3D()
    {
        float area3D = 0.0f;
        for (auto c : charts) area3D += c.second->meshArea;
        return area3D;
    }

    float AreaUV()
    {
        float areaUV = 0.0f;
        for (auto c : charts) areaUV += c.second->ChartArea();
        return areaUV;
    }

    float BorderUV(float *meshBorderLengthUV = nullptr, float *seamLengthUV = nullptr)
    {
        float bl = 0.0f, sl = 0.0f;
        for (auto c : charts) {
            bl += c.second->uvMeshBorder;
            sl += c.second->chartBorder;
        }
        if (meshBorderLengthUV) *meshBorderLengthUV = bl;
        if (seamLengthUV) *seamLengthUV = sl;
        return bl + sl;
    }
};

// Save a copy of the original texture coordinates (this will be used to render the new texture)
// If a face has zero UV space area, sets its tex coord to zero
// TODO is this reasonable? maybe the texture coordinates still convey info about what part of the texture
// data could be relevant
void PreprocessMesh(Mesh &m)
{
    Mesh::PerFaceAttributeHandle<WedgeTexCoordStorage> WTCSh
            = tri::Allocator<Mesh>::GetPerFaceAttribute<WedgeTexCoordStorage>(m, "WedgeTexCoordStorage");
    for (auto &f : m.face) {
        for (int i = 0; i < 3; ++i) {
            if (std::abs(DistortionWedge::AreaUV(&f)) == 0.0f) f.WT(i).P() = Point2f(0.0f, 0.0f);
            WTCSh[&f].tc[i].P() = f.WT(i).P();
            WTCSh[&f].tc[i].N() = f.WT(i).N();
        }
    }
}

// Computes per face connected component ids. Two face with the same id belong
// to the same connected component. Assumes the topology is already updated.
template <class MeshType>
size_t ComputePerFaceConnectedComponentIdAttribute(MeshType &m)
{
    typename MeshType::template PerFaceAttributeHandle<std::size_t> CCIDh
            = tri::Allocator<MeshType>::template GetPerFaceAttribute<std::size_t>(m, "ConnectedComponentID");

    tri::UpdateFlags<MeshType>::FaceClearV(m);
    tri::UpdateTopology<MeshType>::FaceFaceFromTexCoord(m);

    std::stack<typename MeshType::FacePointer> s;
    size_t regionCounter = 0;

    for (auto &f : m.face) {
        if (!f.IsV()) {
            f.SetV();
            s.push(&f);
            size_t id = regionCounter++;
            while (!s.empty()) {
                auto fp = s.top();
                s.pop();
                CCIDh[fp] = id;
                for (int i = 0; i < fp->VN(); ++i) {
                    if (!face::IsBorder(*fp, i)) {
                        auto adj = fp->FFp(i);
                        if (!adj->IsV()) {
                            adj->SetV();
                            s.push(adj);
                        }
                    }
                }

            }
        }
    }
    return regionCounter;
}

template <class MeshType>
std::shared_ptr<ParameterizationData> ComputeParameterizationGraph(
        MeshType &m, std::vector<std::shared_ptr<QImage>> imgVec, float *uvMeshBorder = nullptr)
{
    std::size_t numRegions = ComputePerFaceConnectedComponentIdAttribute<MeshType>(m);

    std::shared_ptr<ParameterizationData> paramData = std::make_shared<ParameterizationData>();
    paramData->textures = imgVec;
    paramData->charts.reserve(numRegions);
    auto CCIDh = tri::Allocator<MeshType>::template GetPerFaceAttribute<std::size_t>(m, "ConnectedComponentID");

    // build parameterization graph
    tri::UpdateTopology<Mesh>::FaceFace(m);
    for (auto &f : m.face) {
        std::size_t regionId = CCIDh[&f];
        paramData->GetChart(regionId)->AddFace(&f, CCIDh);
        // TODO this may be refactored into AddFace
        for (int i = 0; i < f.VN(); ++i) {
            std::size_t adjId = CCIDh[f.FFp(i)];
            if (regionId != adjId) {
                (paramData->GetChart(regionId)->adj).insert(paramData->GetChart(adjId));
            }
        }
    }

    // compute uv mesh border if required
    if (uvMeshBorder) {
        *uvMeshBorder = 0.0f;
        for (auto &f : m.face) {
            for (int i = 0; i < f.VN(); ++i) {
                if (face::IsBorder(f, i)) {
                   *uvMeshBorder += (f.cWT((i+1)%f.VN()).P() - f.cWT(i).P()).Norm();
                }
            }
        }
    }

    return paramData;
}


// Given in input a parameterization graph and two ids, it joins them (by extending id1 with
// the geometry associated to id2)
// ASSUMES FF topology is computed
void MergeRegions(Mesh &m, ParameterizationData &pdata, std::size_t id1, std::size_t id2)
{
    if (id1 == id2) {
        std::cout << "Merging " << id1 << " with itself, nothing to do" << std::endl;
        return;
    }

    //FIXME just select faces using charts[id2]->fpVec
    for (auto fp : pdata.charts[id2]->fpVec) {
        fp->SetS();
    }
/*
    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(m);
    tri::UpdateFlags<Mesh>::FaceClearS(m);

    // select representative face from the second region and expand selection to the whole region
    pdata.charts[id2]->Fp()->SetS();
    tri::UpdateSelection<Mesh>::FaceConnectedFF(m);

    tri::UpdateTopology<Mesh>::FaceFace(m);
*/
    // Build geometry for the poisson mesh FIXME no need for index, use fpVec
    PMesh pm;
    std::vector<std::size_t> index;
    std::vector<bool> sign;
    for (auto &f : m.face) {
        if (f.IsS()) {
            tri::Allocator<PMesh>::AddFace(
                    pm, Point3f(f.WT(0).U(), f.WT(0).V(), 0), Point3f(f.WT(1).U(), f.WT(1).V(), 0), Point3f(f.WT(2).U(), f.WT(2).V(), 0));
            index.push_back(tri::Index(m, f));
            sign.push_back(DistortionWedge::AreaUV(&f) >= 0.0f);
            //std::cout << "pushing face number " << tri::Index(m,f) << std::endl;
        }
    }

    tri::Clean<PMesh>::RemoveDuplicateVertex(pm);
    tri::Allocator<PMesh>::CompactEveryVector(pm);
    tri::UpdateTopology<PMesh>::FaceFace(pm);
    tri::UpdateSelection<PMesh>::Clear(pm);

    auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<std::size_t>(m, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<std::size_t>(m, CHIDh));

    tri::PoissonSolver<PMesh> solver(pm);

    // TODO proper error handling
    if (!solver.IsFeasible()) {
        std::cout << "PoissonSolver.IsFeasible() == false" << std::endl;
        std::exit(-1);
    }

    // shared border between the two regions. this needs to be subtracted
    // from id1 after the merge
    float sharedBorder = 0.0f;

    // compute shared border length and fix its texture coordinates for the poisson solver
    // also subtract uv border length for 3d border edges in order to compute the change in uv space
    // of the parameterized 3d border
    for (std::size_t i = 0; i < index.size(); ++i) {
        PMesh::FaceType &pf = pm.face[i];
        Mesh::FaceType &f = m.face[index[i]];
        assert(f.IsS());

        for (int k = 0; k < f.VN(); ++k) {
            std::size_t adjId = CHIDh[f.FFp(k)];
            if (adjId == id1) { // texture border, copy wedge texture coordinates from adj face
                int j = f.FFi(k);
                TexCoord2f t0 = f.FFp(k)->cWT((j+1)%3);
                TexCoord2f t1 = f.FFp(k)->cWT(j);

                pf.V0(k)->T() = t0;
                pf.V1(k)->T() = t1;

                sharedBorder += (t1.P() - t0.P()).Norm();

                // select vertices (to be fixed later in the poisson solver)
                pf.V0(k)->SetS();
                pf.V1(k)->SetS();
            }
        }
    }

    // poisson solver
    solver.Init();
    solver.SetSelectedAsFixed();
    assert(solver.SolvePoisson());

    // copy texture coordinates back to the original mesh
    assert(index.size() == pm.face.size());
    for (std::size_t i = 0; i < index.size(); ++i) {
        PMesh::FaceType &pf = pm.face[i];
        Mesh::FaceType &f = m.face[index[i]];

        //assert((tri::Distortion<PMesh,true>::AreaUV(&pf) > 0.0f) == sign[i]);

        for (int k = 0; k < f.VN(); ++k) {
            f.WT(k).P() = pf.WT(k).P(); // don't use vertex texcoords, as they are changed after solving the system...
        }
    }

    // Update parameterization data
    auto pc1 = pdata.GetChart(id1);
    auto pc2 = pdata.GetChart(id2);

    pc1->chartBorder -= sharedBorder;

    // clear selection and update id
    for (auto fp : pc2->fpVec) {
        fp->ClearS();
        CHIDh[fp] = pc1->id;
    }

    // readd faces and update adjacency relationships
    for (std::size_t i = 0; i < index.size(); ++i) {
        pc1->AddFace(&m.face[index[i]], CHIDh);
    }

    pc1->adj.erase(pc2);

    for (auto p : pc2->adj) {
        if (p->id != id1) {
            p->adj.erase(pc2);
            p->adj.insert(pc1);
            pc1->adj.insert(p);
        }
    }

    pdata.charts.erase(pc2->id);
}


void ReduceTextureFragmentation(Mesh &m, ParameterizationData &pdata, int maxRegionSize)
{
    using CVal = decltype(pdata.charts)::value_type;
    using BVal = std::unordered_map<std::size_t, float>::value_type;

    assert(maxRegionSize > 0 && maxRegionSize < m.FN());

    tri::UpdateTopology<Mesh>::FaceFace(m);

    while (true) {

        auto t0 = std::chrono::high_resolution_clock::now();

        /* select next region to reduce (rr) POLICY: region with min num of faces
        auto chartEntry = std::min_element(pdata.charts.begin(), pdata.charts.end(),
                                           [](CVal v1, CVal v2) {return v1.second->chartArea > 0.0f && v1.second->FN() < v2.second->FN();});
        assert(chartEntry != pdata.charts.end());

        // if the selected region has FN > threshold*m.FN() STOP
        if (chartEntry->second->FN() > maxRegionSize) {
            std::cout << "Minimal region face count is " << chartEntry->second->FN() << ", stopping." << std::endl;
            break;
        }*/

        // select next island
        auto chartEntry = pdata.charts.begin();
        while(chartEntry != pdata.charts.end()) {
            if (chartEntry->second->adj.size() == 1) break;
            chartEntry++;
        }
        if (chartEntry == pdata.charts.end()) break;

        // select base region as the one adjacent to chartEntry with the longest common border
        std::unordered_map<std::size_t, float> borderMap;
        auto CHIDh  = tri::Allocator<Mesh>::FindPerFaceAttribute<std::size_t>(m, "ConnectedComponentID");
        for (auto fptr : chartEntry->second->fpVec) {
            for (int i = 0; i < fptr->VN(); ++i) {
                auto ffpi = fptr->FFp(i);
                auto adjId = CHIDh[ffpi];
                if (adjId != CHIDh[fptr]) { // texture seam between two regions
                    borderMap[adjId] += DistortionWedge::EdgeLenghtUV(fptr, i);
                }
            }
        }

        auto maxCommonBorder = std::max_element(borderMap.begin(), borderMap.end(), [](BVal v1, BVal v2) {return v1.second < v2.second;});
        assert(maxCommonBorder != borderMap.end());

        // merge the two regions
        std::cout << maxCommonBorder->first << " <- " << chartEntry->first
                  << " (" << chartEntry->second->fpVec.size() << ")" << std::endl;
        MergeRegions(m, pdata, maxCommonBorder->first, chartEntry->first);

        auto t1 = std::chrono::high_resolution_clock::now();

        std::cout << "Iteration took " << std::chrono::duration<float>(t1 - t0).count() << " seconds" << endl;
    }
}

void PrintParameterizationData(std::shared_ptr<ParameterizationData> pdata)
{
    std::cout << pdata->charts.size() << " " << pdata->Area3D() << " "
              << pdata->AreaUV() << " " << pdata->BorderUV() << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " model minFaceCount(int)" << std::endl;
        return -1;
    }

    int minFaceCount = atoi(argv[2]);
    assert(minFaceCount > 0);

    Mesh m;
    int loadmask = 0;

    char *cwd = getcwd(NULL, 0);
    std::string fileName = GetFileName(argv[1]);
    std::string fileDir = GetDirName(argv[1]);
    if (!fileDir.empty()) {
        // change cwd to load materials correctly
        int err = chdir(fileDir.c_str());
        if (err) {
            perror("chdir");
            return -1;
        }
    }

    int r = tri::io::Importer<Mesh>::Open(m, fileName.c_str(), loadmask);
    if (r) {
        std::cerr << tri::io::Importer<Mesh>::ErrorMsg(r) << std::endl;
        return -1;
    }

    if (minFaceCount > m.FN()) {
        std::cout << "minFaceCount > m.FN()" << std::endl;
        return -1;
    }

    std::cout << fileName << " (VN " <<  m.VN() << ", FN " << m.FN() << ")" << std::endl;

    std::vector<std::shared_ptr<QImage>> imgVec;
    for (const string& textureName : m.textures) {
        auto imgptr = std::make_shared<QImage>(textureName.c_str());
        if (imgptr->isNull()) {
            std::cerr << "Unable to load texture file " << textureName.c_str() << std::endl;
            return -1;
        }
        imgVec.push_back(imgptr);
    }

    chdir(cwd);
    cwd = nullptr;

    PreprocessMesh(m);

    tri::RequirePerFaceWedgeTexCoord<Mesh>(m);
    assert(loadmask & tri::io::Mask::IOM_WEDGTEXCOORD);

    float uvMeshBorder;
    
    auto pdata = ComputeParameterizationGraph(m, imgVec, &uvMeshBorder);

    // Print original info
    PrintParameterizationData(pdata);

    auto t0 = std::chrono::high_resolution_clock::now();

    ReduceTextureFragmentation(m, *pdata, minFaceCount);

    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Processing took " << std::chrono::duration<float>(t1 - t0).count() << " seconds" << std::endl;

    std::shared_ptr<QImage> img = RenderTexture(m, imgVec);
    img->save(m.textures[0].c_str(), 0, 100);

    auto pdata2 = ComputeParameterizationGraph(m, imgVec, &uvMeshBorder);
    // Print optimized info
    PrintParameterizationData(pdata2);

    tri::io::ExporterOBJ<Mesh>::Save(m, "out.obj", tri::io::Mask::IOM_WEDGTEXCOORD);

    return 0;
}


