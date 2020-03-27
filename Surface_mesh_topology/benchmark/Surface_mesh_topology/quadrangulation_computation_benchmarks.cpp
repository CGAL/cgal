#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
using namespace CGAL::Surface_mesh_topology;

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  std::string file;
  for (unsigned i=1; i<argc; ++i)
  {
    file = argv[i];
    std::ifstream in(file);
    if (!in.is_open())
    {
      std::cout<<"ERROR reading file " + file<<std::endl;
      exit(EXIT_FAILURE);
    }
    LCC_3_cmap map;
    CGAL::load_off(map, file.c_str());
    Curves_on_surface_topology<LCC_3_cmap> cst(map, true);
    std::cout << "Mesh " << file << std::endl;
    cst.compute_minimal_quadrangulation(true);

    std::cout<<"Initial map: ";
    map.display_characteristics(std::cout)<<std::endl;
    std::cout<<"Reduced map: ";
    cst.get_minimal_quadrangulation().display_characteristics(std::cout)
      <<std::endl;
  }
  return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
