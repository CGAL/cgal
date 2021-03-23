#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_2.h>
#include <CGAL/Geographical_coordinates_traits_2.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Geographical_coordinates_traits_2<K>          Traits;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Traits>    DToS2;

typedef Traits::Point_3                                     Point_3;
typedef Traits::Point_on_sphere_2                           Point;

int main(int argc, char** argv)
{
  std::cout.precision(17);

  const char* filename = (argc > 1) ? argv[1] : "data/poste_france.xyz";

  Traits traits(Point_3(0, 0, 0), 100);


  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Invalid input file: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  typename Traits::Construct_point_on_sphere_2 cps2 = traits.construct_point_on_sphere_2_object();

  std::vector<Point> points;
  double x, y, z;

  while(in >> x >> y >> z)
  {
    Point_3 cp(x, y, z);
    Point ps = cps2(cp);
    std::cout << "Cartesian point: " << cp << " Coordinates on the sphere: " << ps << std::endl;
    points.push_back(ps);
  }

  std::cout << points.size() << " points in input" << std::endl;

  DToS2 dtos(points.begin(), points.end(), traits);

  std::cout << dtos.number_of_vertices() << " vertices" << std::endl;
  std::cout << dtos.number_of_faces() << " solid faces" << std::endl;
  std::cout << dtos.number_of_ghost_faces() << " ghost faces" << std::endl;

  CGAL::write_OFF("result.off", dtos, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
