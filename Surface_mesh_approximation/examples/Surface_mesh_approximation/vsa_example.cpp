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
  // create and read Polyhedron
  Polyhedron input;
  std::ifstream file("data/bear.off");
  if (!file || !(file >> input) || input.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // output data
  Polyhedron output;
  std::vector<std::size_t> triangles;
  std::vector<Kernel::Point_3> anchors;

  // free function interface with named parameters, separated with dots
  CGAL::vsa_mesh_approximation(input, output,
    CGAL::VSA::parameters::init_by_number(200). // seeding by target number of proxies
      init_method(CGAL::VSA_seeding::Hierarchical). // hierarchical init
      iterations(30). // number of relaxation iterations after seeding
      anchor_point(std::back_inserter(anchors)). // get anchor points
      indexed_triangles(std::back_inserter(triangles))); // get indexed triangles

  std::cout << "#anchors: " << anchors.size() << std::endl;
  std::cout << "#triangles: " << triangles.size() << std::endl;

  return EXIT_SUCCESS;
}
