//! \file examples/Arrangement_on_surface_2/conics.cpp
// Constructing an arrangement of various conic arcs.

#include <CGAL/config.h>

#ifdef CGAL_USE_CORE

#include <CGAL/basic.h>
#include <CGAL/draw_arrangement_2.h>

#include "arr_conics.h"
#include "arr_print.h"

int main() {
  Traits traits;
  Arrangement arr(&traits);

  auto ctr_cv = traits.construct_curve_2_object();

  // x-major
  // insert the parabola y = x^2; (-1,1)--(1,1)
  CGAL::insert(arr, ctr_cv(1, 0, 0, 0, -1, 0, CGAL::COUNTERCLOCKWISE,
                           Point(-1, 1), Point(1, 1)));

  // translated
  // Insert the parabola y = x^2 - 2x + 2; (1,1)--(2,2)
  CGAL::insert(arr, ctr_cv(1, 0, 0, -2, -1, 2, CGAL::COUNTERCLOCKWISE,
                           Point(1, 1), Point(2, 2)));
  CGAL::insert(arr, ctr_cv(1, 0, 0, 2, -1, 2, CGAL::COUNTERCLOCKWISE,
                           Point(-2, 2), Point(-1, 1)));

  // rotated
  // Insert the parabola y = x^2 rotated clockwise about theta, such that
  // sin(theta) = 0.6, cos(theta) = 0.8
  CGAL::insert(arr, ctr_cv(16, 9, -24, -15, -20, 0, CGAL::COUNTERCLOCKWISE,
                           Point(Rational(-2,10), Rational(14,10)),
                           Point(Rational(14,10), Rational(2,10))));
  CGAL::insert(arr, ctr_cv(16, 9, 24, 15, -20, 0, CGAL::CLOCKWISE,
                           Point(Rational(2,10), Rational(14,10)),
                           Point(Rational(-14,10), Rational(2,10))));
  CGAL::insert(arr, ctr_cv(16, 9, 24, -15, 20, 0, CGAL::COUNTERCLOCKWISE,
                           Point(Rational(14,10), Rational(-2,10)),
                           Point(Rational(-2,10), Rational(-14,10))));
  CGAL::insert(arr, ctr_cv(16, 9, -24, 15, 20, 0, CGAL::COUNTERCLOCKWISE,
                           Point(Rational(2,10), Rational(-14,10)),
                           Point(Rational(-14,10), Rational(-2,10))));

  // 16*x*x+9*y*y-24*x*y-15*x-20*y

  CGAL::insert(arr, ctr_cv(9, 16, -24, -20, -15, 0, CGAL::COUNTERCLOCKWISE,
                           Point(Rational(2,10), Rational(14,10)),
                           Point(Rational(14,10), Rational(-2,10))));
  CGAL::insert(arr, ctr_cv(9, 16, 24, -20, 15, 0, CGAL::CLOCKWISE,
                           Point(Rational(2,10), Rational(-14,10)),
                           Point(Rational(14,10), Rational(2,10))));
  CGAL::insert(arr, ctr_cv(9, 16, 24, 20, -15, 0, CGAL::COUNTERCLOCKWISE,
                           Point(Rational(-14,10), Rational(-2,10)),
                           Point(Rational(-2,10), Rational(14,10))));
  CGAL::insert(arr, ctr_cv(9, 16, -24, 20, 15, 0, CGAL::COUNTERCLOCKWISE,
                           Point(Rational(-2,10), Rational(-14,10)),
                           Point(Rational(-14,10), Rational(2,10))));

  // 9*x*x+16*y*y-24*x*y+20*x+15*y

  // rotated & translated
  CGAL::insert(arr, ctr_cv(16, 9, -24, -23, -14, 36, CGAL::COUNTERCLOCKWISE,
                           Point(Rational(8,10), Rational(24,10)),
                           Point(Rational(24,10), Rational(12,10))));
  CGAL::insert(arr, ctr_cv(16, 9, 24, 23, -14, 36, CGAL::CLOCKWISE,
                           Point(Rational(-8,10), Rational(24,10)),
                           Point(Rational(-24,10), Rational(12,10))));

  // 16*x*x+9*y*y-24*x*y-23*x-14*y+36

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
