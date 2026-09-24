#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/IO/XYZ.h>

#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

typedef CGAL::Point_set_3<Point> Point_set;

int main ()
{

  Point_set points_only, points_with_normal;

  CGAL::IO::read_XYZ("points_only.xyz", points_only);

  assert(! points_only.has_normal_map());

  CGAL::IO::read_XYZ("points_with_normal.xyz", points_with_normal);

  assert(points_with_normal.has_normal_map());


  return EXIT_SUCCESS;
}
