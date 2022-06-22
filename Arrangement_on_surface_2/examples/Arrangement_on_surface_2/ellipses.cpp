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
  Point p0(0, 0);

  // Insert a full x-major ellipse
  insert(arr, ctr_cv(1, 4, 0, 0, 0, -8));

  // Insert a full y-major ellipse
  insert(arr, ctr_cv(4, 1, 0, 0, 0, -8));

  // Insert a full ellipse (C2), which is (x/4)^2 + (y/2)^2 = 1 rotated by
  // phi = 36.87 degrees (such that sin(phi) = 0.6, cos(phi) = 0.8),
  // yielding: 58x^2 + 72y^2 - 48xy - 360 = 0.
  insert(arr, ctr_cv(58, 72, -48, 0, 0, -360));

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
