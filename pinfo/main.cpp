#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <string>
#include <vector>
#include <memory>
#include <unordered_set>
#include <iostream>
#include <iomanip>

#include <unistd.h>

#include <QImage>

using namespace vcg;

class MeshVertex;
class MeshFace;

struct MeshUsedTypes : public UsedTypes<Use<MeshVertex>::AsVertexType, Use<MeshFace>::AsFaceType> {};

class MeshVertex : public Vertex<MeshUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags> {};
class MeshFace : public Face<MeshUsedTypes, face::VertexRef, face::FFAdj, face::WedgeTexCoord2f, face::BitFlags> {};
class Mesh : public tri::TriMesh<std::vector<MeshVertex>, std::vector<MeshFace>> {};

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
    std::size_t id;
    Mesh::FacePointer fp; // points to a face mapped in the chart

    int faceCount;
    float meshArea;
    float chartArea;

    std::unordered_set<std::shared_ptr<ChartData>> adj;

    ChartData(std::size_t id_, Mesh::FacePointer fp_) : id{id_}, fp{fp_}, faceCount{0}, meshArea{0}, chartArea{0}, adj{} {}

    void Accumulate(const Mesh::FacePointer fptr) {
        faceCount++;
        meshArea += tri::Distortion<Mesh,true>::Area3D(fptr);
        chartArea += tri::Distortion<Mesh,true>::AreaUV(fptr);
    }
};

struct PatameterizationData {
    std::vector<std::shared_ptr<ChartData>> charts;
    std::vector<std::shared_ptr<QImage>> textures;

    PatameterizationData() : charts{}, textures{} {}

    std::shared_ptr<ChartData> CreateChart(const Mesh::FacePointer fptr)
    {
        std::size_t id = charts.size();
        std::shared_ptr<ChartData> chart = std::make_shared<ChartData>(id, fptr);
        charts.push_back(chart);
        return chart;
    }

    std::shared_ptr<ChartData> GetChart(std::size_t i)
    {
        assert(i < charts.size());
        return charts[i];
    }

    float GetArea3D()
    {
        float area3D = 0.0f;
        for (auto c : charts) area3D += c->meshArea;
        return area3D;
    }

    float GetAreaUV()
    {
        float areaUV = 0.0f;
        for (auto c : charts) areaUV += c->chartArea;
        return areaUV;
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

    std::string fileName = GetFileName(argv[1]);
    std::string fileDir = GetDirName(argv[1]);
    if (!fileDir.empty()) {
        // change cwd to load textures correctly...
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

    PatameterizationData paramData;
    std::vector<std::shared_ptr<QImage>> textures;
    for (const string& textureName : m.textures) {
        auto imgptr = std::make_shared<QImage>(textureName.c_str());
        if (imgptr->isNull()) {
            std::cerr << "Unable to load texture file " << textureName.c_str() << std::endl;
            return -1;
        }
        paramData.textures.push_back(imgptr);
    }

    tri::RequirePerFaceWedgeTexCoord<Mesh>(m);
    assert(loadmask & tri::io::Mask::IOM_WEDGTEXCOORD);

    std::cout << "Model info" << std::endl;
    std::cout << "==========" << std::endl;
    std::cout << "Mesh file: " << fileName << " (VN " << m.VN() << ", FN " << m.FN() << ")" << std::endl;
    std::cout << "Texture files:" << std::endl;
    for(std::string &s : m.textures) std::cout << "  " + s << std::endl;

    // Identify each patch (store id as face attribute)

    Mesh::PerFaceAttributeHandle<std::size_t> chartId = tri::Allocator<Mesh>::GetPerFaceAttribute<std::size_t>(m, "ChartID");
    std::stack<Mesh::FacePointer> s;

    tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(m);
    tri::UpdateFlags<Mesh>::FaceClearV(m);

    for (auto fi = m.face.begin(); fi != m.face.end(); ++fi) {
        if (!fi->IsV()) {
            fi->SetV();
            s.push(&*fi);
            std::shared_ptr<ChartData> cdata = paramData.CreateChart(&*fi);
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
        std::size_t chId = chartId[fi];
        for (int i = 0; i < 3; ++i) {
            std::size_t adjId = chartId[fi->FFp(i)];
            if (chId != adjId) {
                (paramData.GetChart(chId)->adj).insert(paramData.GetChart(adjId));
            }
        }
    }

    float tot3D = paramData.GetArea3D();
    float totUV = paramData.GetAreaUV();

    std::cout << std::endl;
    std::cout << "Parameterization info" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "Number of charts: " << paramData.charts.size() << std::endl;
    std::cout << std::setw(6) << "ID" << " " << std::setw(8) << "faces" << " "
              << std::setw(12) << "perc_3D" << " " << std::setw(12) << "perc_UV" << " " << std::setw(12) << "resolution" << " "
              << std::setw(8) << "n_adj" << std::endl;

    for (auto ci = paramData.charts.begin(); ci != paramData.charts.end(); ++ci) {
        std::shared_ptr<QImage> img;
        float texelArea = -0.0f;

        int tid = (*ci)->fp->WT(0).N(); // TODO assumes the whole chart is mapped to the same texture image
        if (tid >= 0) {
            img = paramData.textures[tid];
            texelArea = 1.0f / (img->width() * img->height());
        }

        std::cout << std::setw(6) << (*ci)->id << " " << std::setw(8) << (*ci)->faceCount << " " << std::setw(12) << (*ci)->meshArea / tot3D << " "
                  << std::setw(12) << (*ci)->chartArea / totUV << " " << std::setw(12);
        if (img) std::cout << (*ci)->chartArea / texelArea;
        else std::cout << "*";

        std::cout << " " << std::setw(8) << ((*ci)->adj).size() << std::endl;
    }

    return 0;
}
