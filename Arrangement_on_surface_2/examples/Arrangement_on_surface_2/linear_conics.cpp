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
  Point p0(0, 0);
  Point p1(-1, 0);
  Point p2(0, -1);
  Point p3(0, 1);
  Point p4(1, 0);
  Point p5(Rational(1,2),Rational(1,2));
  Point p6(Rational(-1,2),Rational(1,2));
  Rat_point rp0(0, 0);
  Rat_point rp1(1, 0);
  Rat_point rp2(0, 1);
  Rat_point rp3(0, -1);

  // horizontal
  // insert the segment (0, 0)--(1, 0).
  CGAL::insert(arr, ctr_cv(Rat_segment(rp0, rp1)));
  // insert the segment (0, 0)--(-1, 0).
  CGAL::insert(arr, ctr_cv(0, 0, 0, 0, 1, 0, CGAL::COLLINEAR, p0, p1));

  // vertical
  // insert the segment (0, -1)--(0, 0).
  CGAL::insert(arr, ctr_cv(Rat_segment(rp3, rp0)));

  // translated
  // insert the segment (0, -1)--(1, 0).
  CGAL::insert(arr, ctr_cv(Rat_segment(rp3, rp1)));
  // insert the segment (0, -1)--(-1, 0).
  CGAL::insert(arr, ctr_cv(0, 0, 0, -1, -1, -1, CGAL::COLLINEAR, p2, p1));

  // Special segments
  // horizontal special segment
  CGAL::insert(arr, ctr_cv(p5, p6));

  // vertical special segment
  CGAL::insert(arr, ctr_cv(p0, p3));

  // special translated
  CGAL::insert(arr, ctr_cv(p1, p3));
  CGAL::insert(arr, ctr_cv(p3, p4));

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
