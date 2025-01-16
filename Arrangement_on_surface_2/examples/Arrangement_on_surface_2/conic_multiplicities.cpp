//! \file examples/Arrangement_on_surface_2/conic_multiplicities.cpp
// Handling intersection points with multiplicity between conic arcs.

#include <CGAL/config.h>

#ifdef CGAL_USE_CORE

#include <CGAL/basic.h>
#include <CGAL/Arr_naive_point_location.h>

#include "arr_conics.h"
#include "arr_print.h"

using Naive_pl = CGAL::Arr_naive_point_location<Arrangement>;

int main() {
  Traits traits;
  auto ctr_cv = traits.construct_curve_2_object();

  Arrangement arr(&traits);
  Naive_pl pl(arr);

  // Insert a hyperbolic arc, supported by the hyperbola y = x^2/(1-x)
  // (or: x^2 + xy - y = 0) with the endpoints (-1, 1/2) and (1/2, 1/2).
  // Note that the arc is counterclockwise oriented.
  Point ps1(-1, Rational(1,2));
  Point pt1(Rational(1,2), Rational(1,2));
  Conic_arc cv1 = ctr_cv(1, 0, 1, 0, -1, 0, CGAL::COUNTERCLOCKWISE, ps1, pt1);
  // insert(arr, cv1, pl);

#if 0
  // Insert the bottom half of the circle centered at (0, 1/2) whose radius
  // is 1/2 (therefore its squared radius is 1/4).
  Rat_circle circ2(Rat_point(0, Rational(1,2)), Rational(1,4));
  Point ps2(-Rational(1,2), Rational(1,2));
  Point pt2(Rational(1,2), Rational(1,2));
  Conic_arc cv2 = ctr_cv(circ2, CGAL::COUNTERCLOCKWISE, ps2, pt2);
  insert(arr, cv2, pl);
#endif

  print_arrangement(arr);
  return 0;
}

#else

#include <iostream>

int main ()
{
  std::cout << "Sorry, this example needs GMP and CORE\n";
  return 0;
}

#endif
