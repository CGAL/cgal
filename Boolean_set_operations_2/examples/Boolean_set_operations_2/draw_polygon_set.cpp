/*! \file draw_polygon_set.cpp
 * Drawing a polygon set.
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/draw_polygon_set_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                                Polygon_2;
typedef CGAL::Polygon_with_holes_2<K>                     Polygon_with_holes_2;
typedef CGAL::Polygon_set_2<K>                            Polygon_set_2;
typedef CGAL::Point_2<K>                                  Point_2;

Polygon_2 rectangle(int l)
{
  // Create a rectangle with given side length.
  Polygon_2 P;
  P.push_back(Point_2(-l,-l));
  P.push_back(Point_2(l,-l));
  P.push_back(Point_2(l,l));
  P.push_back(Point_2(-l,l));

  return P;
}

int main()
{
  // Create a large rectangle A, with a hole and a smaller rectangle
  // B inside A's hole.
  Polygon_with_holes_2 A(rectangle(3));
  Polygon_2 H(rectangle(2));
  H.reverse_orientation();
  A.add_hole(H);
  Polygon_2 B(rectangle(1));

  // Add them to a polygon set and draw it.
  Polygon_set_2 S;
  S.insert(A);
  S.insert(B);

  CGAL::draw(S);

  return 0;
}
