// Constructing an arrangement of polycurves.

#include <CGAL/config.h>

#ifndef CGAL_USE_CORE

#include <iostream>

int main() {
  std::cout << "Sorry, this example needs CORE ...\n";
  return 0;
}

#else

#include <vector>

#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_polycurve_traits_2.h>

#include "arr_circular.h"
#include "arr_print.h"

using Polycurve_traits = CGAL::Arr_polycurve_traits_2<Traits>;
using X_monotone_polycurve = Polycurve_traits::X_monotone_curve_2;
using Polycurve = Polycurve_traits::Curve_2;
using Circle_2 = Kernel::Circle_2;
using Polycurve_circ_arc_arrangment = CGAL::Arrangement_2<Polycurve_traits>;

int main() {
  Polycurve_traits traits;
  auto ctr_curve = traits.construct_curve_2_object();
  auto ctr_xcurve = traits.construct_x_monotone_curve_2_object();

  //Containers to store conic curves that will be used to create polycurve.
  std::vector<Curve> curves;
  std::vector<X_monotone_curve> x_curves;

  // Create a circular arc of the circle, directed clockwise from
  // (-1, 0) to (1, 0) centered at (0,0). Note that we orient the
  // supporting circle accordingly.
  Rational_point c1(0, 0);
  Point s1(Number_type(-1), Number_type(0));
  Point t1(Number_type(1), Number_type(0));
  Curve circ_arc1(c1, 1, CGAL::CLOCKWISE, s1, t1);
  curves.push_back(circ_arc1);

  // Create a circular arc of the unit circle, directed clockwise from
  // (1, 0) to (3, 0) centered at (3,0). Note that we orient the
  // supporting circle accordingly.
  Rational_point c2(3, 0);
  Point s2(Number_type(1), Number_type(0));
  Point t2(Number_type(5), Number_type(0));
  Curve circ_arc2(c2, 2, CGAL::CLOCKWISE, s2, t2);
  curves.push_back(circ_arc2);

  // Create polycurve
  Polycurve polycurve_1 = ctr_curve(curves.begin(), curves.end());

  // Empty the vector in order to create another polycurve.
  curves.clear();

  // Create a circular arc of the circle, directed clockwise from
  // (-10, 13) to (-7, 10) centered at (-10,10). Note that we orient the
  // supporting circle accordingly.
  Rational_point c3(-10, 10);
  Point s3(Number_type(-10), Number_type(13));
  Point t3(Number_type(-7), Number_type(10));
  Curve circ_arc3(c3, 3, CGAL::CLOCKWISE, s3, t3);
  curves.push_back(circ_arc3);

  Rational_point c4(-20, 10);
  Point s4(Number_type(-7), Number_type(10));
  Point t4(Number_type(-20), Number_type(23));
  Curve circ_arc4(c4, 13, CGAL::CLOCKWISE, s4, t4);
  curves.push_back(circ_arc4);

  Rational_point c5(-20, 25);
  Point s5(Number_type(-20), Number_type(23));
  Point t5(Number_type(-20), Number_type(27));
  Curve circ_arc5(c5, 2, CGAL::CLOCKWISE, s5, t5);
  curves.push_back(circ_arc5);

  Polycurve polycurve_2 = ctr_curve(curves.begin(), curves.end());

  //circle to be used by x-monotone polycurve
  Rational_point circle_center1(Number_type(10), Number_type(10));
  Circle_2 circ_1(circle_center1, 4, CGAL::CLOCKWISE);
  Point s6(Number_type(8), Number_type(10));
  Point t6(Number_type(12), Number_type(10));
  X_monotone_curve xc1(circ_1, s6, t6, CGAL::CLOCKWISE);
  x_curves.push_back(xc1);

  Rational_point circle_center2(Number_type(13), Number_type(10));
  Circle_2 circ_2(circle_center2, 1, CGAL::CLOCKWISE);
  Point s7(Number_type(12), Number_type(10));
  Point t7(Number_type(14), Number_type(10));
  X_monotone_curve xc2(circ_2, s7, t7, CGAL::CLOCKWISE);
  x_curves.push_back(xc2);

  //create x-monotone polycurve
  X_monotone_polycurve x_polycurve_1 =
    ctr_xcurve(x_curves.begin(), x_curves.end());

  // Insert polycurves to Arrangement and print.
  Polycurve_circ_arc_arrangment polycurve_arrangment(&traits);
  insert(polycurve_arrangment, polycurve_1);
  insert(polycurve_arrangment, polycurve_2);
  insert(polycurve_arrangment, x_polycurve_1);
  std::cout << "Arrangement Statistics:\n";
  print_arrangement(polycurve_arrangment);

  return 0;
}
#endif
