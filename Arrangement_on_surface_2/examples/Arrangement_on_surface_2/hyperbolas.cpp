//! \file examples/Arrangement_on_surface_2/conics.cpp
// Constructing an arrangement of various conic arcs.

#include <CGAL/config.h>

#ifdef CGAL_USE_CORE

#include <CGAL/draw_arrangement_2.h>

#include "arr_conics.h"
#include "arr_print.h"

int main() {
  Traits traits;
  Arrangement arr(&traits);

  auto ctr_cv = traits.construct_curve_2_object();

  // Insert a hyperbolic arc (C1), supported by the hyperbola y = 1/x
  // (or: xy - 1 = 0) with the endpoints (1/4, 4) and (2, 1/2).
  // The arc is counterclockwise oriented.
  CGAL::insert(arr, ctr_cv(0, 0, 1, 0, 0, -1, CGAL::COUNTERCLOCKWISE,
                           Point(Rational(1,4), 4), Point(2, Rational(1,2))));
  CGAL::insert(arr, ctr_cv(0, 0, -1, 0, 0, -1, CGAL::CLOCKWISE,
                           Point(Rational(-1,4), 4), Point(-2, Rational(1,2))));

  CGAL::insert(arr, ctr_cv(2, -1, 0, 0, 0, -2, CGAL::COUNTERCLOCKWISE,
                           Point(3, 4), Point(1, 0)));
  CGAL::insert(arr, ctr_cv(2, -1, 0, 0, 0, -2, CGAL::COUNTERCLOCKWISE,
                           Point(1, 0), Point(3, -4)));
  CGAL::insert(arr, ctr_cv(2, -1, 0, 0, 0, -2, CGAL::CLOCKWISE,
                           Point(-3, 4), Point(-1, 0)));
  CGAL::insert(arr, ctr_cv(2, -1, 0, 0, 0, -2, CGAL::CLOCKWISE,
                           Point(-1, 0), Point(-3, -4)));

  CGAL::insert(arr, ctr_cv(-1, 2, 0, 0, 0, -2, CGAL::CLOCKWISE,
                           Point(4, 3), Point(0, 1)));
  CGAL::insert(arr, ctr_cv(-1, 2, 0, 0, 0, -2, CGAL::CLOCKWISE,
                           Point(0, 1), Point(-4, 3)));
  CGAL::insert(arr, ctr_cv(-1, 2, 0, 0, 0, -2, CGAL::COUNTERCLOCKWISE,
                           Point(4, -3), Point(0, -1)));
  CGAL::insert(arr, ctr_cv(-1, 2, 0, 0, 0, -2, CGAL::COUNTERCLOCKWISE,
                           Point(0, -1), Point(-4, -3)));

  CGAL::insert(arr, ctr_cv(4, 46, -144, 0, 0, -100, CGAL::COUNTERCLOCKWISE,
                           Point(-5, 0),
                           Point(Rational(14, 10), Rational(48, 10))));
  CGAL::insert(arr, ctr_cv(4, 46, -144, 0, 0, -100, CGAL::COUNTERCLOCKWISE,
                           Point(5, 0),
                           Point(Rational(-14, 10), Rational(-48, 10))));
  // 4*x*x + 46*y*y - 144*x*y - 100

  CGAL::insert(arr, ctr_cv(46, 4,  -144, 0, 0, -100, CGAL::CLOCKWISE,
                           Point(0, -5),
                           Point(Rational(48, 10), Rational(14, 10))));
  CGAL::insert(arr, ctr_cv(46, 4,  -144, 0, 0, -100, CGAL::CLOCKWISE,
                           Point(0, 5),
                           Point(Rational(-48, 10), Rational(-14, 10))));
  // 46*x*x + 4*y*y - 144*x*y - 100

  CGAL::insert(arr, ctr_cv(4, 46, 144, 0, 0, -100, CGAL::CLOCKWISE,
                           Point(-5, 0),
                           Point(Rational(14,10), Rational(-48,10))));
  CGAL::insert(arr, ctr_cv(4, 46, 144, 0, 0, -100, CGAL::CLOCKWISE,
                           Point(5, 0),
                           Point(Rational(-14,10), Rational(48,10))));
  // 4*x*x + 46*y*y + 144*x*y - 100

  CGAL::insert(arr, ctr_cv(46, 4, 144, 0, 0, -100, CGAL::COUNTERCLOCKWISE,
                           Point(0, -5),
                           Point(Rational(-48,10), Rational(14,10))));
  CGAL::insert(arr, ctr_cv(46, 4, 144, 0, 0, -100, CGAL::COUNTERCLOCKWISE,
                           Point(0, 5),
                           Point(Rational(48,10), Rational(-14,10))));

  print_arrangement_size(arr);

  CGAL::draw(arr);

  return 0;
}

#else

#include <iostream>

int main() {
  std::cout << "Sorry, this example needs GMP and CORE\n";
  return 0;
}

#endif
