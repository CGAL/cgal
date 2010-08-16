#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>

#include <iostream>
#include <iterator>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                  Point_2;
typedef CGAL::Random_points_in_square_2<
     Point_2,
     CGAL::Creator_uniform_2< double, Point_2 > >
                                    Point_generator;
int main() {
  // create 500-gon and write it into a window:
  CGAL::random_convex_set_2(
            500,
            std::ostream_iterator<Point_2>(std::cout, "\n"),
            Point_generator( 0.5));
  return 0;
}
