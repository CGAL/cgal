#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

namespace VSA = CGAL::Surface_mesh_approximation;

int main(int argc, char** argv)
{
  const char* filename = (argc > 1) ? argv[1] : "data/bear.off";

  // reads input surface triangle mesh
  Mesh mesh;
  if(!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, mesh) ||
     !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

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
