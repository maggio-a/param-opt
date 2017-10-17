#include <iostream>

#include "dcpsolver.h"

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

using namespace vcg;

int main( int argc, char **argv )
{
    Mesh m;
    if(argc < 2 ) return -1;
    printf("Reading %s\n",argv[1]);
    int loadMask;
    int ret = tri::io::Importer<Mesh>::Open(m,argv[1], loadMask);
    if (ret != 0)
    {
        printf("Unable to open %s for '%s'\n",argv[1],tri::io::Importer<Mesh>::ErrorMsg(ret));
        return -1;
    }
    printf("Mesh has %i vn %i fn\n",m.VN(),m.FN());

    DCPSolver solver{m};
    bool solved = solver.Solve();

    if (!solved) std::cout << "Not solved" << std::endl;

    //tri::UpdateTexture<MyMesh>::WedgeTexFromVertexTex(m);
    tri::io::ExporterPLY<Mesh>::Save(m,"out.ply",tri::io::Mask::IOM_WEDGTEXCOORD);

    return 0;
}
