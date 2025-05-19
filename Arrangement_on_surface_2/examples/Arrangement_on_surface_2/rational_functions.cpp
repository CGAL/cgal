//! \file examples/Arrangement_2/rational_functions.cpp
// Constructing an arrangement of arcs of rational functions.

#include <CGAL/config.h>

#ifndef CGAL_USE_CORE
#include <iostream>

int main() {
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return 0;
}

#else

#include "arr_rat_functions.h"
#include "arr_print.h"

int main() {
  CGAL::IO::set_pretty_mode(std::cout);             // for nice printouts.

  // Define a traits class object and a constructor for rational functions.
  Traits traits;
  auto construct = traits.construct_x_monotone_curve_2_object();

  // Define a polynomial representing x.
  Polynomial x = CGAL::shift(Polynomial(1), 1);

  // Define a container storing all arcs.
  std::vector<Traits::X_monotone_curve_2> arcs;

  // Create an arc (C1) supported by the polynomial y = x^4 - 6x^2 + 8,
  // defined over the (approximate) interval [-2.1, 2.1].
  Polynomial P1 = CGAL::ipower(x,4) - 6*x*x + 8;
  Alg_real l(Bound(-2.1)), r(Bound(2.1));
  arcs.push_back(construct(P1, l, r));

  // Create an arc (C2) supported by the function y = x / (1 + x^2),
  // defined over the interval [-3, 3].
  Polynomial P2 = x;
  Polynomial Q2 = 1 + x*x;
  arcs.push_back(construct(P2, Q2, Alg_real(-3), Alg_real(3)));

  // Create an arc (C3) supported by the parbola y = 8 - x^2,
  // defined over the interval [-2, 3].
  Polynomial P3 = 8 - x*x;
  arcs.push_back(construct(P3, Alg_real(-2), Alg_real(3)));

  // Create an arc (C4) supported by the line y = -2x,
  // defined over the interval [-3, 0].
  Polynomial P4 = -2*x;
  arcs.push_back(construct(P4, Alg_real(-3), Alg_real(0)));

  // Construct the arrangement of the four arcs.
  Arrangement arr(&traits);
  insert(arr, arcs.begin(), arcs.end());
  print_arrangement(arr);

  return 0;
}

#endif
