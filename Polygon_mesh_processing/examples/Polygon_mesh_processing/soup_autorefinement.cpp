#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
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
  std::vector<std::array<std::size_t, 3>> input_triangles;
  CGAL::IO::read_polygon_soup(filename, input_points, input_triangles);

  std::vector<Point> output_points;
  std::vector<std::array<std::size_t, 3>> output_triangles;
  PMP::autorefine_soup_output(input_points, input_triangles, output_points, output_triangles);

  CGAL::IO::write_polygon_soup("autorefined.off", output_points, output_triangles, CGAL::parameters::stream_precision(17));

  return 0;
}
