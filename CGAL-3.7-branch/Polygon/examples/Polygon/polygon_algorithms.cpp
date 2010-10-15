#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
using std::cout; using std::endl;

void check_inside(Point pt, Point *pgn_begin, Point *pgn_end, K traits)
{
  cout << "The point " << pt;
  switch(CGAL::bounded_side_2(pgn_begin, pgn_end,pt, traits)) {
    case CGAL::ON_BOUNDED_SIDE :
      cout << " is inside the polygon.\n";
      break;
    case CGAL::ON_BOUNDARY:
      cout << " is on the polygon boundary.\n";
      break;
    case CGAL::ON_UNBOUNDED_SIDE:
      cout << " is outside the polygon.\n";
      break;
  }
}

int main()
{
  Point points[] = { Point(0,0), Point(5.1,0), Point(1,1), Point(0.5,6)};

  // check if the polygon is simple.
  cout << "The polygon is "
    << (CGAL::is_simple_2(points, points+4, K()) ? "" : "not ")
    << "simple." << endl;

  check_inside(Point(0.5, 0.5), points, points+4, K());
  check_inside(Point(1.5, 2.5), points, points+4, K());
  check_inside(Point(2.5, 0), points, points+4, K());

  return 0;
}
