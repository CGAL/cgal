#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/triangle_soup_boolean_operations_3.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename1 =
    (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");
  const std::string filename2 =
    (argc > 2) ? argv[2] : CGAL::data_file_path("meshes/eight.off");

  std::vector<Point> points1, points2;
  std::vector<std::array<std::size_t, 3>> triangles1, triangles2;

  if(!CGAL::IO::read_polygon_soup(filename1, points1, triangles1) ||
     !CGAL::IO::read_polygon_soup(filename2, points2, triangles2))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::vector<Point> points_output;
  std::vector<std::array<std::size_t, 3>> triangles_output;

  PMP::compute_union(points1, triangles1,
                     points2, triangles2,
                     points_output, triangles_output);

  std::cout << "Union was successfully computed\n";

  CGAL::IO::write_polygon_soup(
    "union.off", points_output, triangles_output,
    CGAL::parameters::stream_precision(17));

  return 0;
}