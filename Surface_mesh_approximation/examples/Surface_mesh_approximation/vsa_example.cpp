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
  Polyhedron mesh;
  std::ifstream input("data/bear.off");
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // output data
  Polyhedron out_mesh;
  std::vector<int> tris;
  std::vector<Kernel::Point_3> anchor_pos;

  // free function interface with named parameters
  CGAL::vsa_mesh_approximation(mesh, out_mesh,
    CGAL::VSA::parameters::number_of_proxies(200). // number of fitting proxies
      number_of_iterations(30). // number of iterations
      init_method(1). // hierarchical init
      anchor_point(std::back_inserter(anchor_pos)). // get anchor points
      indexed_triangles(std::back_inserter(tris))); // get indexed triangles

  std::cout << "#anchor_pos " << anchor_pos.size() << std::endl;
  std::cout << "#tris " << tris.size() << std::endl;

  return EXIT_SUCCESS;
}
