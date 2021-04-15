#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

namespace VSA = CGAL::Surface_mesh_approximation;

int main()
{
  // read input surface triangle mesh
  Mesh mesh;
  std::ifstream file("data/bear.off");
  file >> mesh;

  // The output will be an indexed triangle mesh
  std::vector<Kernel::Point_3> anchors;
  std::vector<std::array<std::size_t, 3> > triangles;

  // free function interface with named parameters
  VSA::approximate_triangle_mesh(mesh,
    CGAL::parameters::verbose_level(VSA::MAIN_STEPS).
    max_number_of_proxies(200).
    anchors(std::back_inserter(anchors)). // anchor points
    triangles(std::back_inserter(triangles))); // indexed triangles

  std::cout << "#anchor points: " << anchors.size() << std::endl;
  std::cout << "#triangles: " << triangles.size() << std::endl;

  return EXIT_SUCCESS;
}
