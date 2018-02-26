#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/vsa_mesh_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main()
{
  // read input surface triangle mesh 
  Polyhedron input;
  std::ifstream file("data/mask.off");
  if (!file || !(file >> input) || input.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // output polyhedral surface and indexed face set
  Polyhedron output;
  std::vector<Kernel::Point_3> points;
  std::vector<std::vector<std::size_t> > triangles; // triplets of indices

  // free function interface with named parameters
  bool valid_polyhedron = CGAL::vsa_mesh_approximation(input,
    std::back_inserter(points),
    std::back_inserter(triangles),
    CGAL::Surface_mesh_approximation::parameters::seeding_method(CGAL::Hierarchical). // hierarchical seeding
    max_nb_proxies(200). // seeding with maximum number of proxies
    nb_of_iterations(30). // number of clustering iterations after seeding
    // output to polyhedron
    output_mesh(&output)); // valid when the indexed face set represents a 2-manifold, oriented surface

  std::cout << "#anchor points: " << points.size() << std::endl;
  std::cout << "#triangles: " << triangles.size() << std::endl;

  if (valid_polyhedron)
    std::cout << "oriented, 2-manifold output." << std::endl;

  return EXIT_SUCCESS;
}
