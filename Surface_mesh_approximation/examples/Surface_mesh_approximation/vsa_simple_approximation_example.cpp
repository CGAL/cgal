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
  std::vector<std::size_t> triangles;
  std::vector<Kernel::Point_3> points;

  // free function interface with named parameters
  bool valid_polyhedron = CGAL::VSA::mesh_approximation(input,
	  CGAL::VSA::parameters::init_method(CGAL::VSA::Hierarchical). // hierarchical init
	  refine_until_proxies(200). // refine until target number of proxies
	  iterations(30). // number of relaxation iterations after seeding
	  anchor_points(std::back_inserter(points)). // get anchor points
	  indexed_triangles(std::back_inserter(triangles))); // get indexed triangles

  std::cout << "#anchor points: " << points.size() << std::endl;
  std::cout << "#triangles: " << triangles.size() / 3 << std::endl;

  return EXIT_SUCCESS;
}
