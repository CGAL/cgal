#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/hilbert_sort.h>
#include <vector>
#include <array>
#include <string>
#include <filesystem>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Points = std::vector<Point_3>;
using Face = std::array<std::size_t,3>;
using Polygons = std::vector<Face>;




int main(int argc, char* argv[])
{

  const std::string input = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");
  const std::string output = std::filesystem::path(input).stem().string() + "_sorted.off";
  Points points;
  Polygons polygons;
  CGAL::IO::read_polygon_soup(input, points, polygons);

  CGAL::Polygon_mesh_processing::hilbert_sort_polygon_soup(points, polygons);


  CGAL::IO::write_polygon_soup(output, points, polygons);

  return EXIT_SUCCESS;
}

