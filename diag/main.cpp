#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <vector>
#include <memory>
#include <unordered_set>
#include <iostream>
#include <iomanip>

using namespace vcg;

class MeshVertex;
class MeshFace;

struct MeshUsedTypes : public UsedTypes<Use<MeshVertex>::AsVertexType, Use<MeshFace>::AsFaceType> {};

class MeshVertex : public Vertex<MeshUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags> {};
class MeshFace : public Face<MeshUsedTypes, face::VertexRef, face::FFAdj, face::WedgeTexCoord2f, face::BitFlags> {};
class Mesh : public tri::TriMesh<std::vector<MeshVertex>, std::vector<MeshFace>> {};

struct ChartData {
    std::size_t id;
    Mesh::FacePointer fp;

    int faceCount;
    float meshArea;
    float chartArea;

    std::unordered_set<std::shared_ptr<ChartData>> adj;

    ChartData(std::size_t id_, Mesh::FacePointer fp_) : id{id_}, fp{fp_}, faceCount{0}, meshArea{0}, chartArea{0}, adj{} {}

    void Accumulate(const Mesh::FacePointer f) {
        faceCount++;
        meshArea += tri::Distortion<Mesh,true>::Area3D(f);
        chartArea += tri::Distortion<Mesh,true>::AreaUV(f);
    }
};

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "No model specified" << std::endl;
        return -1;
    }

    Mesh m;
    int loadmask = 0;
    int r = tri::io::Importer<Mesh>::Open(m, argv[1], loadmask);
    if (r) {
        std::cout << tri::io::Importer<Mesh>::ErrorMsg(r) << std::endl;
        return -1;
    }

    tri::RequirePerFaceWedgeTexCoord<Mesh>(m);
    assert(loadmask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD);

    std::cout << argv[1] << " loaded" << std::endl;

    // Identify each patch (store id as face attribute)

    std::vector<std::shared_ptr<ChartData>> charts{};
    Mesh::PerFaceAttributeHandle<std::size_t> chartId = tri::Allocator<Mesh>::GetPerFaceAttribute<std::size_t>(m, "ChartID");
    std::stack<Mesh::FacePointer> s;

    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(m);
    tri::UpdateFlags<Mesh>::FaceClearV(m);

    for (auto fi = m.face.begin(); fi != m.face.end(); ++fi) {
        if (!fi->IsV()) {
            fi->SetV();
            s.push(&*fi);
            charts.push_back(std::make_shared<ChartData>(charts.size(), &*fi));
            auto cdata = charts.back();
            while (!s.empty()) {
                auto fp = s.top();
                s.pop();
                chartId[fp] = cdata->id;

                cdata->Accumulate(fp);

                for (int i = 0; i < 3; ++i) {
                    if (!face::IsBorder(*fp, i)) {
                        auto fadj = fp->FFp(i);
                        if (!fadj->IsV()) {
                            fadj->SetV();
                            s.push(&*fadj);
                        }
                    }
                }
            }
        }
    }

    // Compute adjacency list for each chart

    tri::UpdateTopology<Mesh>::FaceFace(m);

    for (auto fi = m.face.begin(); fi != m.face.end(); ++fi) {
        std::size_t faceId = chartId[fi];
        for (int i = 0; i < 3; ++i) {
            std::size_t adjId = chartId[fi->FFp(i)];
            if (faceId != adjId) {
                (charts[faceId]->adj).insert(charts[adjId]);
            }
        }
    }

    std::cout << "Chart info" << std::endl
              << "Parameterization charts: " << charts.size() << std::endl;

    std::cout << std::setw(6) << "ID" << " " << std::setw(8) << "numfaces" << " " << std::setw(12) << "mesh_area" << " " << std::setw(12) << "uv_area" << " adj" << std::endl;
    for (auto ci = charts.begin(); ci != charts.end(); ++ci) {
        std::cout << std::setw(6) << (*ci)->id << " " << std::setw(8) << (*ci)->faceCount << " " << std::setw(12) << (*ci)->meshArea << " "
                  << std::setw(12) << (*ci)->chartArea << " ";
        std::cout << "[ ";
        for (auto neighbor : (*ci)->adj)
            std::cout << neighbor->id << " ";
        std::cout << "]" << std::endl;
    }

    return 0;
}
