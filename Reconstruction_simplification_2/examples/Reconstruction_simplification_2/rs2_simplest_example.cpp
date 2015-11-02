// Example for Reconstruction_simplification_2, with no mass
// attributes for the input points

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Reconstruction_simplification_2.h>

#include <fstream>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;

typedef CGAL::Reconstruction_simplification_2<K>            Rs_2;

int main ()
{
  // Generate a set of random points.
  std::vector<Point> points;
  CGAL::Random_points_on_square_2<Point> point_generator(1.);
  CGAL::cpp11::copy_n(point_generator, 100, std::back_inserter(points));

  Rs_2 rs2(points);
  rs2.run(100); // 100 steps

  return 0;
}
