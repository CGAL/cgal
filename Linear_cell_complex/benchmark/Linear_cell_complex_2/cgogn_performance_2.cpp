//== INCLUDES =================================================================
#include "cgogn_performance_2.h"
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
    CGoGN_performance_2().run(argv[i], "output_cgogn_2.off");
  }
  return 0;
}
//=============================================================================
