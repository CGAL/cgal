#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/vsa_mesh_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main()
{
  // create polyhedral surface and read input surface triangle mesh 
  Polyhedron input;
  std::ifstream file("data/mask.off");
  if (!file || !(file >> input) || input.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // output polyhedral surface and indexed triangle mesh
  Polyhedron output;
  std::vector<std::vector<std::size_t> > triangles; // contains triplets of indices
  std::vector<Kernel::Point_3> points;

  // free function interface with named parameters
  bool valid_polyhedron = CGAL::VSA::mesh_approximation(input,
    std::back_inserter(points),
    std::back_inserter(triangles),
    CGAL::VSA::parameters::init_method(CGAL::VSA::Hierarchical). // hierarchical init
    refine_until_proxies(200). // refine until target number of proxies
    iterations(30). // number of relaxation iterations after seeding
    // output to polyhedron
    output_mesh(output)); // valid only if the indexed triangles construct a 2-manifold and oriented surface

  std::cout << "#anchor points: " << points.size() << std::endl;
  std::cout << "#triangles: " << triangles.size() << std::endl;

  if (valid_polyhedron)
    std::cout << "oriented 2-manifold output." << std::endl;

  return EXIT_SUCCESS;
}
