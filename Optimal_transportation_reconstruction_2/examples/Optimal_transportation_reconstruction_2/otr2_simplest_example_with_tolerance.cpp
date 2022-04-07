// Simplest example with tolerance for
// Optimal_transportation_reconstruction_2, with no mass attributes
// for the input points

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>

#include <fstream>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;

typedef CGAL::Optimal_transportation_reconstruction_2<K>    Otr;

int main ()
{
  // Generate 100 random points on the boundary of a square.
  std::vector<Point> points;
  CGAL::Random_points_on_square_2<Point> point_generator(1.);
  std::copy_n(point_generator, 100, std::back_inserter(points));

  Otr otr(points); // no mass given, one unit mass per point assumed
  otr.run_under_wasserstein_tolerance(0.1);


  return 0;
}
