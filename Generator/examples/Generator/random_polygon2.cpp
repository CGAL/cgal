#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                 Point_2;
typedef std::list<Point_2>                         Container;
typedef CGAL::Polygon_2<K, Container>              Polygon_2;
typedef CGAL::Random_points_in_square_2< Point_2 > Point_generator;

int main() {
  Polygon_2 polygon;
  // create 50-gon and write it into a window:
  CGAL::random_polygon_2(50, std::back_inserter(polygon),
                         Point_generator(0.5));
  std::cout << polygon;
  return 0;
}
