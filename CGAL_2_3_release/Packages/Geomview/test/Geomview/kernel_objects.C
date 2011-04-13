// Simple test program for Geomview_stream.
// See the demo directory for more extensive use.
//
// Sylvain Pion, 2000.

#include <CGAL/Cartesian.h>

#include <CGAL/IO/Geomview_stream.h>

typedef CGAL::Cartesian<double> K;

void test_parse_point()
{
  const char *test_point="( 123 456 789 1 )";
  K::Point_3 p;
  CGAL::parse_point(test_point, p);
  CGAL_assertion(p == K::Point_3(123, 456, 789));
}

int main()
{
  test_parse_point();
  return 0;
}
