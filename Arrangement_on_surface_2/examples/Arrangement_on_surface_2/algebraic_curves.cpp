//! \file examples/Arrangement_on_surface_2/algebraic_curves.cpp
// Constructing an arrangement of algebraic curves.

#include <iostream>
#include <CGAL/config.h>

#if (!CGAL_USE_CORE) && (!CGAL_USE_LEDA) && (!(CGAL_USE_GMP && CGAL_USE_MPFI))
int main ()
{
  std::cout << "Sorry, this example needs CORE, LEDA, or GMP+MPFI ..."
            << std::endl;
  return 0;
}

#else

#include <CGAL/basic.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>

#include "integer_type.h"
#include "arr_print.h"

using Traits = CGAL::Arr_algebraic_segment_traits_2<Integer>;
using Arrangement = CGAL::Arrangement_2<Traits>;
using Curve = Traits::Curve_2;
using Polynomial = Traits::Polynomial_2;

int main() {
  CGAL::IO::set_pretty_mode(std::cout);             // for nice printouts.

  Traits traits;
  Arrangement arr(&traits);

  // Functor to create a curve from a Polynomial.
  auto construct_cv = traits.construct_curve_2_object();

  Polynomial x = CGAL::shift(Polynomial(1), 1, 0);
  Polynomial y = CGAL::shift(Polynomial(1), 1, 1);

  // Construct an unbounded line (C1) with the equation 3x - 5y - 2 = 0.
  Polynomial f1 = 3*x - 5*y - 2;
  Curve cv1 = construct_cv(f1);
  std::cout << "Inserting curve " << f1 << std::endl;
  CGAL::insert(arr, cv1);

  // Construct the ellipse (C2) with the equation x^2 + 3*y^2 - 10 = 0.
  Polynomial f2 = CGAL::ipower(x, 2) + 3*CGAL::ipower(y, 2) - 10;
  Curve cv2 = construct_cv(f2);
  std::cout << "Inserting curve " << f2 << std::endl;
  CGAL::insert(arr, cv2);

  // Construct a cubic curve (C3) with the equation x^2 + y^2 + xy^2 = 0,
  // with isolated point at (0,0) and vertical asymptote at x = 1.
  Polynomial f3 = CGAL::ipower(x, 2) + CGAL::ipower(y, 2) +
    x*CGAL::ipower(y, 2);
  Curve cv3 = construct_cv(f3);
  std::cout << "Inserting curve " << f3 << std::endl;
  CGAL::insert(arr, cv3);

  // Construct a curve of degree 6 (C4) with the equation
  // x^6 + y^6 - x^3y^3 - 12 = 0.
  Polynomial f4 = CGAL::ipower(x, 6) + CGAL::ipower(y, 6) -
    CGAL::ipower(x, 3)*CGAL::ipower(y, 3) - 12;
  Curve cv4 = construct_cv(f4);
  std::cout << "Inserting curve " << f4 << std::endl;
  CGAL::insert(arr, cv4);

  print_arrangement_size(arr);                // print the arrangement size
  return 0;
}

#endif
