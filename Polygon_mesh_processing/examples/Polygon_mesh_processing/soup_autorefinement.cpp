#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point;

typedef CGAL::Surface_mesh<Point>                           Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const std::string filename = argc == 1 ? CGAL::data_file_path("meshes/elephant.off")
                                         : std::string(argv[1]);

  std::vector<Point> input_points;
  std::vector<std::vector<std::size_t> > input_triangles;

  CGAL::IO::read_polygon_soup(filename, input_points, input_triangles);
  PMP::repair_polygon_soup(input_points, input_triangles);

  Mesh mesh;
  PMP::orient_polygon_soup(input_points, input_triangles);
  PMP::polygon_soup_to_polygon_mesh(input_points, input_triangles, mesh);

  PMP::autorefine(mesh);

  CGAL::IO::write_polygon_mesh("autorefined.off", mesh, CGAL::parameters::stream_precision(17));

  return 0;
}
