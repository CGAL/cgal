#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/helpers.h>

#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  for(int i = 1; i < argc; ++i) {
    Surface_mesh mesh;
    if(!PMP::IO::read_polygon_mesh(argv[i], mesh)) {
      std::cerr << "Invalid input: " << argv[i] << '\n';
    }
    if(!is_triangle_mesh(mesh)) {
      std::cout << argv[i] << std::endl;
    }
  }

  return 0;
}
