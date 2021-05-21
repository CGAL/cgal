#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_2.h>
#include <CGAL/Projection_on_sphere_traits_3.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Projection_on_sphere_traits_3<K>              Traits;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Traits>    DToS2;

typedef DToS2::Point_3                                      Point;

int main(int argc, char** argv)
{
  std::cout.precision(17);

  const char* filename = (argc > 1) ? argv[1] : "data/radar.xyz";

  std::vector<Point> points;
  double x, y, z;

  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Invalid input file: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  while(in >> x >> y >> z)
    points.emplace_back(x, y, z);

  std::cout << points.size() << " points in input" << std::endl;

  Traits traits(Point(0, 0, 0), 100);
  DToS2 dtos(points.begin(), points.end(), traits);

  std::cout << dtos.number_of_vertices() << " vertices" << std::endl;
  std::cout << dtos.number_of_faces() << " solid faces" << std::endl;

  CGAL::IO::write_OFF("result.off", dtos, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
