#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_approximation.h>

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

  // output indexed triangles
  std::vector<Kernel::Point_3> points;
  std::vector<CGAL::cpp11::array<std::size_t, 3> > triangles; // triplets of indices

  // free function interface with named parameters
  bool is_manifold = CGAL::mesh_approximation(input,
    CGAL::Surface_mesh_approximation::parameters::seeding_method(CGAL::Hierarchical). // hierarchical seeding
    max_nb_proxies(200). // seeding with maximum number of proxies
    nb_of_iterations(30). // number of clustering iterations after seeding
    anchor_points(std::back_inserter(points)). // anchor points
    indexed_triangles(std::back_inserter(triangles))); // indexed triangles

  std::cout << "#anchor points: " << points.size() << std::endl;
  std::cout << "#triangles: " << triangles.size() << std::endl;

  if (is_manifold)
    std::cout << "oriented, 2-manifold output." << std::endl;

  return EXIT_SUCCESS;
}
