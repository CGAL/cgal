//== INCLUDES =================================================================
#include "cgogn_performance_3.h"
#include "lcc_performance_3.h"
#include "openvolumemesh_performance.h"
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
      std::cout << "LCC_3\t" << std::endl;
      LCC_performance_3().run(argv[i], "output_lcc_3.tetmesh");
    }

    {
      std::cout << "CGoGN_3\t" << std::endl;
      CGoGN_performance_3().run(argv[i], "output_cgogn_3.tetmesh");
    }

    {
      std::string filename(argv[i]);
      size_t pos = filename.rfind(".");    // position of "." in filename
      filename.replace (pos, std::string::npos, ".ovm");
      std::cout << "OpenVolumeMesh\t" << std::endl;
      OpenVolumeMesh_performance().run(filename.c_str(), "output_openvolumemesh.ovm");
    }
  }

  return 0;
}
//=============================================================================
