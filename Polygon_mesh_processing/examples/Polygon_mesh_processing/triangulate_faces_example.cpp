#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                    Point;
typedef CGAL::Surface_mesh<Point>          Surface_mesh;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/cube_quad.off";
  std::ifstream input(filename);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  CGAL::Polygon_mesh_processing::triangulate_faces(mesh);

  CGAL_assertion(is_pure_triangle(mesh));

  std::ofstream cube_off("cube_tri.off");
  cube_off << mesh;
  cube_off.close();

  return 0;
}
