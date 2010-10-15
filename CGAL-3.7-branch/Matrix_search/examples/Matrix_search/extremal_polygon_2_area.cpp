#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/extremal_polygon_2.h>
#include <iostream>
#include <vector>

typedef double                                    FT;

typedef CGAL::Cartesian<FT>                       Kernel;

typedef Kernel::Point_2                           Point;
typedef std::vector<int>                          Index_cont;
typedef CGAL::Polygon_2<Kernel>                   Polygon_2;
typedef CGAL::Random_points_in_square_2<Point>    Generator;

int main() {

  int n = 10;
  int k = 5;

  // generate random convex polygon:
  Polygon_2 p;
  CGAL::random_convex_set_2(n, std::back_inserter(p), Generator(1));
  std::cout << "Generated Polygon:\n" << p << std::endl;

  // compute maximum area incribed k-gon of p:
  Polygon_2 k_gon;
  CGAL::maximum_area_inscribed_k_gon_2(
    p.vertices_begin(), p.vertices_end(), k, std::back_inserter(k_gon));
  std::cout << "Maximum area " << k << "-gon:\n"
            << k_gon << std::endl;

  return 0;
}
