/*! \file draw_polygon_with_holes_2.cpp
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>

using K=CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2=K::Point_2;
using Polygon_2=CGAL::Polygon_2<K>;
using Polygon_with_holes_2=CGAL::Polygon_with_holes_2<K>;

Polygon_2 triangle(double x, double y)
{ // Create a triangle (x, y)
  Polygon_2 P;
  P.push_back(Point_2(x, y));
  P.push_back(Point_2(x+0.5, y+1));
  P.push_back(Point_2(x+1, y));
  return P;
}

Polygon_2 rectangle(double l)
{   // Create a rectangle with given side length.
  Polygon_2 P;
  P.push_back(Point_2(0, 0));
  P.push_back(Point_2(l, 0));
  P.push_back(Point_2(l, l));
  P.push_back(Point_2(0, l));
  return P;
}

int main()
{
  // Create a large rectangle A, with 3 triangle holes
  Polygon_with_holes_2 A(rectangle(4));
  Polygon_2 H1(triangle(1, 1));
  Polygon_2 H2(triangle(2, 1));
  Polygon_2 H3(triangle(1.5, 2));
  A.add_hole(H1);
  A.add_hole(H2);
  A.add_hole(H3);
  CGAL::draw(A);

  return 0;
}
