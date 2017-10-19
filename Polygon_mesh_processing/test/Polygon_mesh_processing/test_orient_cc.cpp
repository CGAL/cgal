#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <iostream>
#include <fstream>

using namespace CGAL;
namespace PMP = Polygon_mesh_processing;

typedef Exact_predicates_inexact_constructions_kernel Kernel;
typedef Surface_mesh<Kernel::Point_3> Surface_mesh;


int main()
{

  std::ifstream input("data-coref/nested_cubes_invalid_volume.off");
  assert(input);
  ::Surface_mesh sm;
  input >> sm;
  PMP::orient_connected_components(sm, PMP::parameters::all_default());
  return 0;
}
