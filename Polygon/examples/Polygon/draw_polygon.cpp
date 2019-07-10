#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_polygon_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef CGAL::Point_2<K>                                    Point;

int main()
{
  // create a polygon and put some points in it
  Polygon_2 p;
  p.push_back(Point(0,0));
  p.push_back(Point(4,0));
  p.push_back(Point(4,4));
  p.push_back(Point(2,2));
  p.push_back(Point(0,4));

  CGAL::draw(p);

  return EXIT_SUCCESS;
}
