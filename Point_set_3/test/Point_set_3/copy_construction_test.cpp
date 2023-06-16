#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/grid_simplify_point_set.h>

#include <boost/unordered_set.hpp>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;
typedef std::array<unsigned char, 3> Color;

int main (int, char**)
{
  Point_set original;
  original.insert({1, 2, 3});
  original.insert({4, 5, 6});

  Point_set copy_constructed{original};
  assert(copy_constructed.number_of_points() == original.number_of_points());
  assert(copy_constructed.point(0) == original.point(0));
  assert(copy_constructed.point(1) == original.point(1));

  Point_set copy_assigned;
  copy_assigned = original;
  assert(copy_assigned.number_of_points() == original.number_of_points());
  assert(copy_assigned.point(0) == original.point(0));
  assert(copy_assigned.point(1) == original.point(1));

}
