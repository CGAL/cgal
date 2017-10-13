#include <fstream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Set_movable_separability_2/Single_mold_translational_casting/pullout_directions.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;

namespace SMS = CGAL::Set_movable_separability_2;
namespace casting = SMS::Single_mold_translational_casting;

// The main program:
int main(int  argc, char* argv[])
{
  Polygon_2 polygon;

  const char* filename = (argc > 1) ? argv[1] : "polygon.dat";
  std::ifstream input_file(filename);
  if (! input_file.is_open()) {
    std::cerr << "Failed to open the " << filename << std::endl;
    return -1;
  }
  input_file >> polygon;
  input_file.close();

  // Example for pullout_directions_single_mold_translational_casting_2
  size_t index(0);
  for (auto e_it = polygon.edges_begin(); e_it != polygon.edges_end(); ++e_it,
         ++index)
  {
    auto res = casting::pullout_directions(polygon, e_it);
    if (res.first) {
      std::cout << "The polygon is castable using edge " << index
                << " in range " << res.second.first
                << " to " << res.second.second << std::endl;
    }
    else {
      std::cout << "The polygon is not castable using edge " << index
                << std::endl;
    }
  }

  return 0;
}
