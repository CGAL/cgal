//== INCLUDES =================================================================
#include "polyhedron_performance.h"
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
    Polyhedron_performance().run(argv[i], "output_polyhedron.off");
  }

  return 0;
}
//=============================================================================
