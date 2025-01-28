//! \file examples/Arrangement_on_surface_2/circular_arc.cpp
// Constructing an arrangement of various circular arcs and line segments.

#include <CGAL/draw_arrangement_2.h>

#include "arr_circular.h"
#include "arr_print.h"

int main() {
  std::list<Curve> curves;

  // Create a circle (C1) centered at the origin with squared radius 2.
  curves.push_back(Curve(Circle(Rational_point(0, 0), Number_type(2))));

  // Create a circle (C2) centered at (2, 3) with radius 3/2. Note that
  // as the radius is rational we use a different curve constructor.
  Number_type three_halves = Number_type(3) / Number_type(2);
  curves.push_back(Curve(Rational_point(2, 3), three_halves));

  // Create a segment (C3) of the line (y = x) with rational endpoints.
  Segment s3 = Segment(Rational_point(-2, -2), Rational_point(2, 2));
  curves.push_back(Curve(s3));

  // Create a line segment (C4) with the same supporting line (y = x), but
  // having one endpoint with irrational coordinates.
  CoordNT sqrt_15 = CoordNT(0, 1, 15); // = sqrt(15)
  curves.push_back(Curve(s3.supporting_line(),
                         Point(3, 3), Point(sqrt_15, sqrt_15)));

  // Create a circular arc (C5) that is the upper half of the circle centered at
  // (1, 1) with squared radius 3. Create the circle with clockwise orientation,
  // so the arc is directed from (1 - sqrt(3), 1) to (1 + sqrt(3), 1).
  Rational_point c5(1, 1);
  Circle circ5(c5, 3, CGAL::CLOCKWISE);
  CoordNT one_minus_sqrt_3(1, -1, 3);
  CoordNT one_plus_sqrt_3(1, 1, 3);
  Point s5(one_minus_sqrt_3, CoordNT(1));
  Point t5(one_plus_sqrt_3, CoordNT(1));
  curves.push_back(Curve(circ5, s5, t5));

  // Create an arc (C6) of the unit circle, directed clockwise from
  // (-1/2, sqrt(3)/2) to (1/2, sqrt(3)/2).
  // The supporting circle is oriented accordingly.
  Rational_point c6(0, 0);
  Number_type half = Number_type(1) / Number_type(2);
  CoordNT sqrt_3_div_2(Number_type(0), half, 3);
  Point s6(-half, sqrt_3_div_2);
  Point t6(half, sqrt_3_div_2);
  curves.push_back(Curve(c6, 1, CGAL::CLOCKWISE, s6, t6));

  // Create a circular arc (C7) defined by two endpoints and a midpoint,
  // all having rational coordinates. This arc is the upper right
  // quarter of a circle centered at the origin with radius 5.
  curves.push_back(Curve(Rational_point(0, 5), Rational_point(3, 4),
                         Rational_point(5, 0)));

  // Construct the arrangement of the curves and print its size.
  Arrangement  arr;
  insert(arr, curves.begin(), curves.end());
  print_arrangement_size(arr);
  CGAL::draw(arr);
  return 0;
}
