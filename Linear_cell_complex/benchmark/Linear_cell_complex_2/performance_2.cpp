//== INCLUDES =================================================================
#include "polyhedron_performance.h"
#include "lcc_performance_2.h"
#include "openmesh_performance.h"
#include "surface_mesh_performance.h"
#include "cgogn_performance_2.h"
//=============================================================================
int main(int argc, char** argv)
{
  if (argc < 2)
    {
        std::cerr << "Usage:\nperformance <input-mesh>\n";
        exit(1);
    }

  for (int i=1; i<argc; ++i)
  {
    std::cout<<"**************** "<<argv[i]<<" ****************"<<std::endl;
    {
      std::cout << "Polyhedron\t" << std::endl;
      Polyhedron_performance().run(argv[i], "output_polyhedron.off");
    }

    {
      std::cout << "LCC_2\t" << std::endl;
      LCC_performance_2().run(argv[1], "output_lcc_2.off");
    }

    {
      std::cout << "OpenMesh\t" << std::endl;
      OpenMesh_performance().run(argv[1], "output_open_mesh.off");
    }

    {
      std::cout << "Surface_mesh\t" << std::endl;
      Surface_mesh_performance().run(argv[1], "output_surface_mesh.off");
    }

    {
      std::cout << "CGoGn_2\t" << std::endl;
      CGoGN_performance_2().run(argv[1], "output_cgogn_2.off");
    }
  }

  return 0;
}
//=============================================================================
