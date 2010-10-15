//! \file examples/Arrangement_on_surface_2/conic_multiplicities.cpp
// Handling intersection points with multiplicity between conic arcs.

#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl; 
  return 0;
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>

#include "arr_print.h"

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef Rat_kernel::Point_2                           Rat_point_2;
typedef Rat_kernel::Segment_2                         Rat_segment_2;
typedef Rat_kernel::Circle_2                          Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, 
                                 Alg_kernel,
                                 Nt_traits>           Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Conic_arc_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2> Naive_pl;

int main ()
{
  Arrangement_2  arr;
  Naive_pl       pl (arr);

  // Insert a hyperbolic arc, supported by the hyperbola y = x^2/(1-x)
  // (or: x^2 + xy - y = 0) with the endpoints (-1, 1/2) and (1/2, 1/2).
  // Note that the arc is counterclockwise oriented.
  Point_2        ps1 (-1, Rational(1,2));
  Point_2        pt1 (Rational(1,2), Rational(1,2));
  Conic_arc_2    cv1 (1, 0, 1, 0, -1, 0, CGAL::COUNTERCLOCKWISE, ps1, pt1);

  insert (arr, cv1, pl);

  // Insert the bottom half of the circle centered at (0, 1/2) whose radius
  // is 1/2 (therefore its squared radius is 1/4).
  Rat_circle_2   circ2 (Rat_point_2(0, Rational(1,2)), Rational(1,4));
  Point_2        ps2 (-Rational(1,2), Rational(1,2));
  Point_2        pt2 (Rational(1,2), Rational(1,2));
  Conic_arc_2    cv2 (circ2, CGAL::COUNTERCLOCKWISE, ps2, pt2);
  
  insert (arr, cv2, pl);

  // Print the resulting arrangement.
  print_arrangement (arr);

  return 0;
}

#endif
