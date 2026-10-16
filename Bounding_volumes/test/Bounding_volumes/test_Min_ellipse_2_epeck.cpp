#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <cassert>
#include <vector>

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = K::Point_2;

void check(const std::vector<Point>& points, int support)
{
  CGAL::Min_ellipse_2<CGAL::Min_ellipse_2_traits_2<K>> ellipse(points.begin(), points.end());
  assert(ellipse.is_valid());
  assert(ellipse.number_of_support_points() == support);
  for(const Point& p : points)
    assert(!ellipse.has_on_unbounded_side(p));
  ellipse.insert(Point(0, 0));
  assert(ellipse.is_valid());
}

int main()
{
  check({Point(0, 0), Point(1, 0), Point(0, 1)}, 3);
  check({Point(-1, -1), Point(-1, 1), Point(1, 1), Point(1, -1)}, 4);
  check({Point(5, 0), Point(3, 4), Point(-3, 4), Point(-5, 0), Point(0, -5)}, 5);
  check({Point(0, 0), Point(1, 1), Point(2, 2)}, 2);
}
