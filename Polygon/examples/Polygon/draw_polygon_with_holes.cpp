#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_with_holes_2<K>                       Polygon_with_holes_2;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef CGAL::Point_2<K>                                    Point;

int main()
{
  // create a polygon with three holes
  Polygon_2 outer_polygon;
  outer_polygon.push_back(Point(0,0)); outer_polygon.push_back(Point(9,0));
  outer_polygon.push_back(Point(6,8)); outer_polygon.push_back(Point(5,3));
  outer_polygon.push_back(Point(2,8)); outer_polygon.push_back(Point(0,8));

  std::vector<Polygon_2> holes(3);
  holes[0].push_back(Point(6,2)); holes[0].push_back(Point(7,1));
  holes[0].push_back(Point(7,3)); holes[0].push_back(Point(6,3));
  holes[0].push_back(Point(5,2));

  holes[1].push_back(Point(2,1)); holes[1].push_back(Point(3,1));
  holes[1].push_back(Point(3,3)); holes[1].push_back(Point(2,2));
  holes[1].push_back(Point(1,2));

  holes[2].push_back(Point(1,4)); holes[2].push_back(Point(2,4));
  holes[2].push_back(Point(2,5)); holes[2].push_back(Point(3,5));
  holes[2].push_back(Point(3,6)); holes[2].push_back(Point(1,6));

  Polygon_with_holes_2 p(outer_polygon, holes.begin(), holes.end());

  // And draw it.
  CGAL::draw(p);

  return EXIT_SUCCESS;
}
