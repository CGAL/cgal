//! \file examples/Arrangement_2/unbounded_rational_functions.cpp
// Constructing an arrangement of unbounded portions of rational functions.

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
  AK1 ak1;
  Traits traits(&ak1);
  auto construct = traits.construct_curve_2_object();

  // Define a polynomial representing x.
  Polynomial x = CGAL::shift(Polynomial(1), 1);

  // Define a container storing all arcs.
  std::vector<Traits::Curve_2> arcs;

  // Create the arcs (C1, C'1) of the rational functions (y = 1 / x, y = -1 / x).
  Polynomial P1(1);
  Polynomial minusP1(-P1);
  Polynomial Q1 = x;
  arcs.push_back(construct(P1, Q1));
  arcs.push_back(construct(minusP1, Q1));

  // Create the bounded segments (C2, C'2) of the parabolas (y = -4*x^2 + 3)
  // and (y = 4*x^2 - 3), defined over [-sqrt(3)/2, sqrt(3)/2].
  Polynomial P2 = -4*x*x+3;
  Polynomial minusP2 = -P2;
  std::vector<std::pair<Alg_real, int> > roots;
  ak1.solve_1_object()(P2, std::back_inserter(roots));// [-sqrt(3)/2, sqrt(3)/2]
  arcs.push_back(construct(P2, roots[0].first, roots[1].first));
  arcs.push_back(construct(minusP2, roots[0].first, roots[1].first));

  // Create the arcs (C3, C'3) of (i) the rational function (y = 1 / 2*x) for
  // x > 0, and (ii) the rational function (y = -1 / 2*x) for x < 0.
  Polynomial P3(1);
  Polynomial minusP3(-P3);
  Polynomial Q3 = 2*x;
  arcs.push_back(construct(P3, Q3, Alg_real(0), true));
  arcs.push_back(construct(minusP3, Q3, Alg_real(0), false));

  // Construct the arrangement of the six arcs and print its size.
  Arrangement arr(&traits);
  insert(arr, arcs.begin(), arcs.end());
  print_unbounded_arrangement_size(arr);

  return 0;
}

#endif
