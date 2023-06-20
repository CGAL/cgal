#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/point_generators_3.h>

#include <boost/unordered_set.hpp>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;

int main() {

  Point_set points;
  assert(points.number_of_points() == 0);

  // Add a large number of points using a generator
  std::size_t reserved_points_count = 1000;
  CGAL::Random_points_in_cube_3<Point> generator;
  points.reserve(reserved_points_count);
  assert(points.number_of_points() == 0);
  for (std::size_t i = 0; i < reserved_points_count; ++i)
    points.insert(*(generator++));
  assert(points.number_of_points() == reserved_points_count);

  // Add more points without making a reservation beforehand
  std::size_t additional_points_count = 100;
  for (std::size_t i = 0; i < additional_points_count; ++i)
    points.insert(*(generator++));
  assert(points.number_of_points() == reserved_points_count + additional_points_count);

}
