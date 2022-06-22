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

  // Insert a hyperbolic arc (C1), supported by the hyperbola y = 1/x
  // (or: xy - 1 = 0) with the endpoints (1/4, 4) and (2, 1/2).
  // The arc is counterclockwise oriented.
  insert(arr, ctr_cv(0, 0, 1, 0, 0, -1, CGAL::COUNTERCLOCKWISE,
                     Point(Rational(1,4), 4), Point(2, Rational(1,2))));

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
