#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/vsa_mesh_segmentation.h>

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

  // TODO: get output indexed face set

  // free function interface with named parameters
  CGAL::VSA::mesh_segmentation(input,
    CGAL::VSA::parameters::init_method(CGAL::VSA::Hierarchical). // hierarchical init
    max_nb_proxies(200). // refine until target number of proxies
    iterations(30)); // number of relaxation iterations after seeding

  return EXIT_SUCCESS;
}
