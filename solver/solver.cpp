#include <iostream>

#include "dcpsolver.h"

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

using namespace vcg;

int main( int argc, char **argv )
{
    PMesh m;
    if(argc < 2 ) return -1;
    printf("Reading %s\n",argv[1]);
    int loadMask;
    int ret = tri::io::Importer<PMesh>::Open(m,argv[1], loadMask);
    if (ret != 0)
    {
        printf("Unable to open %s for '%s'\n",argv[1],tri::io::Importer<PMesh>::ErrorMsg(ret));
        return -1;
    }
    printf("Mesh has %i vn %i fn\n",m.VN(),m.FN());

    DCPSolver<PMesh> solver{m};
    bool solved = solver.Solve();

    if (!solved) std::cout << "Not solved" << std::endl;

    int n1 = 0;
    int n2 = 0;
    int n3 = 0;

    for (auto &f : m.face) {
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

    tri::io::ExporterPLY<PMesh>::Save(m,"out.ply",tri::io::Mask::IOM_WEDGTEXCOORD);

    return 0;
}
