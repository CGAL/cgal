#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <iostream>

struct Kernel : public CGAL::Cartesian<double> {};

typedef Kernel::Point_2                           Point_2;
typedef Kernel::Line_2                            Line_2;
typedef CGAL::Polygon_2<Kernel>                   Polygon_2;
typedef CGAL::Random_points_in_square_2<Point_2>  Generator;

int main()
{
  // build a random convex 20-gon p
  Polygon_2 p;
  CGAL::random_convex_set_2(20, std::back_inserter(p), Generator(1.0));
  std::cout << p << std::endl;

  // compute the minimal enclosing strip p_m of p
  Line_2 p_m[2];
  CGAL::min_strip_2(p.vertices_begin(), p.vertices_end(), p_m);
  std::cout << p_m[0] << "\n" << p_m[1] << std::endl;

  return 0;
}
