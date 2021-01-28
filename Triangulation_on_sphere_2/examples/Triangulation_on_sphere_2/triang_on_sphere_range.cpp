#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_2.h>
#include <CGAL/Projection_sphere_traits_3.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Projection_sphere_traits_3<K>                 Traits;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Traits>    DToS2;

typedef Traits::Point_3                                     Point;

int main(int argc, char** argv)
{
  std::cout.precision(17);

  const char* filename = (argc > 1) ? argv[1] : "poste_france.data";

  int n;
  std::vector<Point> points;
  double x, y, z;

  std::ifstream in(filename);
  in >> n;
  while(in >> x >> y >> z)
    points.emplace_back(x, y, z);

  std::cout << points.size() << " points in input" << std::endl;

  Traits traits(Point(0, 0, 0), 100);
  DToS2 dtos(points.begin(), points.end(), traits);

  std::ofstream out("result.off");
  out.precision(17);
  CGAL::write_OFF(out, dtos);
}
