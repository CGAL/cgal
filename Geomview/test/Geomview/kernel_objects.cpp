// Simple test program for Geomview_stream.
// See the demo directory for more extensive use.
//
// Sylvain Pion, 2000.

#include <CGAL/basic.h>

#ifndef CGAL_USE_GEOMVIEW
#include <iostream>
int main()
{
  std::cout << "Geomview untested on Windows." << std::endl;
  return 0;
}
#else

#include <CGAL/Cartesian.h>

#include <CGAL/IO/Geomview_stream.h>

#include <cassert>

typedef CGAL::Cartesian<double> K;

void test_parse_point()
{
  const char *test_point="( 123 456 789 1 )";
  double x, y, z, w;
  CGAL::Geomview_stream::parse_point(test_point, x, y, z, w);
  K::Point_3 p(x, y, z, w);
  assert(p == K::Point_3(123, 456, 789));
}

int main()
{
  test_parse_point();
  return 0;
}
#endif
