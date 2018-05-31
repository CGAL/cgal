#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/approximate_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main()
{
  // read input surface triangle mesh 
  Polyhedron input;
  std::ifstream file("data/mask.off");
  file >> input;

  // output indexed triangles
  std::vector<Kernel::Point_3> anchors;
  std::vector<CGAL::cpp11::array<std::size_t, 3> > triangles; // triplets of indices

  // free function interface with named parameters
  bool is_manifold = CGAL::approximate_mesh(input,
    CGAL::VSA::parameters::seeding_method(CGAL::Hierarchical). // hierarchical seeding
    max_nb_proxies(200). // seeding with maximum number of proxies
    nb_of_iterations(30). // number of clustering iterations after seeding
    anchors(std::back_inserter(anchors)). // anchor vertices
    triangles(std::back_inserter(triangles))); // indexed triangles

  std::cout << "#anchor vertices: " << anchors.size() << std::endl;
  std::cout << "#triangles: " << triangles.size() << std::endl;

  if (is_manifold)
  {
    std::cout << "oriented, 2-manifold output." << std::endl;
    // TODO: convert from soup to polyhedron mesh
  }

  return EXIT_SUCCESS;
}
