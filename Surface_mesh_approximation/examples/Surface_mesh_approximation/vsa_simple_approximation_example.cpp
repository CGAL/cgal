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
  file >> input;

  // The output will be an indexed triangle mesh
  std::vector<Kernel::Point_3> points;
  std::vector<CGAL::cpp11::array<std::size_t, 3> > triangles;

  // free function interface with named parameters
  CGAL::mesh_approximation(input,
    CGAL::Surface_mesh_approximation::parameters::max_nb_proxies(200).
    anchor_points(std::back_inserter(points)). // anchor points
    indexed_triangles(std::back_inserter(triangles))); // indexed triangles

  std::cout << "#vertices: " << points.size() << std::endl;
  std::cout << "#triangles: " << triangles.size() << std::endl;

  return EXIT_SUCCESS;
}
