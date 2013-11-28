//== INCLUDES =================================================================
#include "surface_mesh_performance.h"
//=============================================================================
int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage:\n" << argv[0] << " <input-mesh>\n";
        exit(1);
    }

  for (int i=1; i<argc; ++i)
  {
    std::cout<<"**************** "<<argv[i]<<" ****************"<<std::endl;
    Surface_mesh_performance().run(argv[1], "output_surface_mesh.off");
  }

  return 0;
}
//=============================================================================
