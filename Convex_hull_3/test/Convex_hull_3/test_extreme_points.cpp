#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <vector>
#include <cassert>

typedef CGAL::Exact_rational                          NT;
typedef CGAL::Cartesian<NT>                           K;
typedef CGAL::Convex_hull_traits_3<K>                 Traits;
typedef Traits::Polygon_mesh                          Polyhedron_3;
typedef K::Point_3                                    Point_3;

void test_function_overload()
{
  std::vector<Point_3> points;
  points.push_back(Point_3(0,0,0));
  points.push_back(Point_3(10,0,0));
  points.push_back(Point_3(0,10,0));
  points.push_back(Point_3(0,0,10));
  points.push_back(Point_3(5,5,5));
  points.push_back(Point_3(2,5,3));
  points.push_back(Point_3(1,3,2));

  Polyhedron_3 polyhedron;
  CGAL::convex_hull_3(points.begin(), points.end(), polyhedron);

  std::vector<Point_3> extreme_points;
  CGAL::extreme_points_3(points, std::back_inserter(extreme_points));

  typedef Polyhedron_3::Point_iterator Point_iterator;
  int i = 0;
  for (Point_iterator p_it = polyhedron.points_begin(); p_it != polyhedron.points_end(); ++p_it)
  {
    Point_3 p = *p_it;
    CGAL_assertion(p == extreme_points[i]);
    i++;
  }
}



int main()
{

  test_function_overload();

  return 0;
}
