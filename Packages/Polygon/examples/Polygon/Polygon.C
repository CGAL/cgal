// file: examples/Polygon/Polygon.C

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <iostream>

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon;
using std::cout; using std::endl;


int main()
{
  Point points[] = { Point(0,0), Point(5.1,0), Point(1,1), Point(0.5,6)};
  Polygon pgn(points, points+4);

  // check if the polygon is simple.
  cout << "The polygon is " << 
    (pgn.is_simple() ? "" : "not ") << "simple." << endl;

  // check if the polygon is convex
  cout << "The polygon is " << 
    (pgn.is_convex() ? "" : "not ") << "convex." << endl;

  return 0;
}

