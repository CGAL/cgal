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

  // Insert a full x-major ellipse
  CGAL::insert(arr, ctr_cv(1, 4, 0, 0, 0, -16, CGAL::COUNTERCLOCKWISE,
                           Point(4,0), Point(0,2)));
  CGAL::insert(arr, ctr_cv(1, 4, 0, 0, 0, -16, CGAL::COUNTERCLOCKWISE,
                           Point(0,2), Point(-4,0)));
  CGAL::insert(arr, ctr_cv(1, 4, 0, 0, 0, -16, CGAL::COUNTERCLOCKWISE,
                           Point(-4,0), Point(0,-2)));
  CGAL::insert(arr, ctr_cv(1, 4, 0, 0, 0, -16, CGAL::COUNTERCLOCKWISE,
                           Point(0,-2), Point(4,0)));

  // Insert a full y-major ellipse
  CGAL::insert(arr, ctr_cv(4, 1, 0, 0, 0, -16));

  // Insert the full ellipse (x/4)^2 + (y/2)^2 = 1 clockwise rotated by
  // phi = 36.87 degrees (such that sin(phi) = 0.6, cos(phi) = 0.8),
  CGAL::insert(arr, ctr_cv(52, 73, 72, 0, 0, -400));

  // Insert the full ellipse (x/4)^2 + (y/2)^2 = 1 counter clockwise rotated by
  // phi = 36.87 degrees (such that sin(phi) = 0.6, cos(phi) = 0.8),
  CGAL::insert(arr, ctr_cv(52, 73, -72, 0, 0, -400));

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
