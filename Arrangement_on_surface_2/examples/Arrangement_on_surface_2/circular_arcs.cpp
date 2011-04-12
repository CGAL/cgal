//! \file examples/Arrangement_on_surface_2/circular_arc.cpp
// Constructing an arrangement of various circular arcs and line segments.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef Kernel::Circle_2                              Circle_2;
typedef Kernel::Segment_2                             Segment_2;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>     Traits_2;
typedef Traits_2::CoordNT                             CoordNT;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Curve_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

int main ()
{
  std::list<Curve_2>  curves;

  // Create a circle centered at the origin with squared radius 2.
  Kernel::Point_2  c1 = Kernel::Point_2 (0, 0);
  Circle_2         circ1 = Circle_2 (c1, Number_type (2));
  
  curves.push_back (Curve_2 (circ1));

  // Create a circle centered at (2,3) with radius 3/2 - note that
  // as the radius is rational we use a different curve constructor.
  Kernel::Point_2  c2 = Kernel::Point_2 (2, 3);
  
  curves.push_back (Curve_2 (c2, Number_type(3, 2)));

  // Create a segment of the line (y = x) with rational endpoints.
  Kernel::Point_2  s3 = Kernel::Point_2 (-2, -2);
  Kernel::Point_2  t3 = Kernel::Point_2 (2, 2);
  Segment_2        seg3 = Segment_2 (s3, t3);

  curves.push_back (Curve_2 (seg3));

  // Create a line segment with the same supporting line (y = x), but
  // having one endpoint with irrational coefficients.
  CoordNT          sqrt_15 = CoordNT (0, 1, 15); // = sqrt(15)
  Point_2          s4 = Point_2 (3, 3);
  Point_2          t4 = Point_2 (sqrt_15, sqrt_15);

  curves.push_back (Curve_2 (seg3.supporting_line(), s4, t4));

  // Create a circular arc that correspond to the upper half of the
  // circle centered at (1,1) with squared radius 3. We create the
  // circle with clockwise orientation, so the arc is directed from
  // (1 - sqrt(3), 1) to (1 + sqrt(3), 1).
  Kernel::Point_2  c5 = Kernel::Point_2 (1, 1);
  Circle_2         circ5 = Circle_2 (c5, 3, CGAL::CLOCKWISE);
  CoordNT          one_minus_sqrt_3 = CoordNT (1, -1, 3);
  CoordNT          one_plus_sqrt_3 = CoordNT (1, 1, 3);
  Point_2          s5 = Point_2 (one_minus_sqrt_3, CoordNT (1));
  Point_2          t5 = Point_2 (one_plus_sqrt_3, CoordNT (1));

  curves.push_back (Curve_2 (circ5, s5, t5));

  // Create a circular arc of the unit circle, directed clockwise from
  // (-1/2, sqrt(3)/2) to (1/2, sqrt(3)/2). Note that we orient the
  // supporting circle accordingly.
  Kernel::Point_2  c6 = Kernel::Point_2 (0, 0);
  CoordNT          sqrt_3_div_2 = CoordNT (Number_type(0), Number_type(1,2), Number_type(3));
  Point_2          s6 = Point_2 (Number_type (-1, 2), sqrt_3_div_2);
  Point_2          t6 = Point_2 (Number_type (1, 2), sqrt_3_div_2);
  
  curves.push_back (Curve_2 (c6, 1, CGAL::CLOCKWISE, s6, t6));

  // Create a circular arc defined by two endpoints and a midpoint,
  // all having rational coordinates. This arc is the upper-right
  // quarter of a circle centered at the origin with radius 5.
  Kernel::Point_2  s7 = Kernel::Point_2 (0, 5);
  Kernel::Point_2  mid7 = Kernel::Point_2 (3, 4);
  Kernel::Point_2  t7 = Kernel::Point_2 (5, 0);

  curves.push_back (Curve_2 (s7, mid7, t7));

  // Construct the arrangement of the curves.
  Arrangement_2    arr;

  insert (arr, curves.begin(), curves.end());
  
  // Print the size of the arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  return 0;
}
