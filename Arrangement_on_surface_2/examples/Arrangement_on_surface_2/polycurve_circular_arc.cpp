// Constructing an arrangement of polycurves.

#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE

#include <iostream>

int main()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return 0;
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <vector>
#include <list>
#include <CGAL/Arr_polycurve_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include "arr_print.h"

typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>     Arc_traits_2;
typedef CGAL::Arr_polycurve_traits_2<Arc_traits_2>    Polycurve_arc_traits_2;

typedef Arc_traits_2::CoordNT                         CoordNT;
typedef Arc_traits_2::Point_2                         Point_2;
typedef Arc_traits_2::Curve_2                         Arc_section_2;
typedef Arc_traits_2::X_monotone_curve_2              Arc_x_monotone_section_2;

typedef Polycurve_arc_traits_2::X_monotone_curve_2    X_monotone_polycurve;
typedef Polycurve_arc_traits_2::Curve_2               Polycurve;
typedef Kernel::Circle_2                              Circle_2;
typedef CGAL::Arrangement_2<Polycurve_arc_traits_2>
  Polycurve_circ_arc_arrangment;

int main()
{
  Polycurve_arc_traits_2 traits;

  //Containers to store conic curves that will be used to create polycurve.
  std::vector<Arc_section_2> curves;
  std::vector<Arc_x_monotone_section_2> x_curves;

  // Create a circular arc of the circle, directed clockwise from
  // (-1, 0) to (1, 0) centered at (0,0). Note that we orient the
  // supporting circle accordingly.
  Kernel::Point_2 c1(0, 0);
  Point_2 s1(Number_type(-1, 1), Number_type(0, 1));
  Point_2 t1(Number_type(1, 1), Number_type(0, 1));
  Arc_section_2 circ_arc1(c1, 1, CGAL::CLOCKWISE, s1, t1);
  curves.push_back(circ_arc1);

  // Create a circular arc of the unit circle, directed clockwise from
  // (1, 0) to (3, 0) centered at (3,0). Note that we orient the
  // supporting circle accordingly.
  Kernel::Point_2 c2(3, 0);
  Point_2 s2(Number_type(1, 1), Number_type(0, 1));
  Point_2 t2(Number_type(5, 1), Number_type(0, 1));
  Arc_section_2 circ_arc2(c2, 2, CGAL::CLOCKWISE, s2, t2);
  curves.push_back(circ_arc2);

  // Create polycurve
  Polycurve polycurve_1 =
    traits.construct_curve_2_object()(curves.begin(), curves.end());

  // Empty the vector in order to create another polycurve.
  curves.clear();

  // Create a circular arc of the circle, directed clockwise from
  // (-10, 13) to (-7, 10) centered at (-10,10). Note that we orient the
  // supporting circle accordingly.
  Kernel::Point_2 c3(-10, 10);
  Point_2 s3(Number_type(-10, 1), Number_type(13, 1));
  Point_2 t3(Number_type(-7, 1), Number_type(10, 1));
  Arc_section_2 circ_arc3(c3, 3, CGAL::CLOCKWISE, s3, t3);
  curves.push_back(circ_arc3);

  Kernel::Point_2 c4(-20, 10);
  Point_2 s4(Number_type(-7, 1), Number_type(10, 1));
  Point_2 t4(Number_type(-20, 1), Number_type(23, 1));
  Arc_section_2 circ_arc4(c4, 13, CGAL::CLOCKWISE, s4, t4);
  curves.push_back(circ_arc4);

  Kernel::Point_2 c5(-20, 25);
  Point_2 s5(Number_type(-20, 1), Number_type(23, 1));
  Point_2 t5(Number_type(-20, 1), Number_type(27, 1));
  Arc_section_2 circ_arc5(c5, 2, CGAL::CLOCKWISE, s5, t5);
  curves.push_back(circ_arc5);

  Polycurve polycurve_2 =
    traits.construct_curve_2_object()(curves.begin(), curves.end());

  //circle to be used by x-monotone polycurve
  Kernel::Point_2 circle_center1(Number_type(10, 1), Number_type(10, 1));
  Circle_2 circ_1(circle_center1, 4, CGAL::CLOCKWISE);
  Point_2 s6(Number_type(8, 1), Number_type(10, 1));
  Point_2 t6(Number_type(12, 1), Number_type(10, 1));
  Arc_x_monotone_section_2 xc1(circ_1, s6, t6, CGAL::CLOCKWISE);
  x_curves.push_back(xc1);

  Kernel::Point_2 circle_center2(Number_type(13, 1), Number_type(10, 1));
  Circle_2 circ_2(circle_center2, 1, CGAL::CLOCKWISE);
  Point_2 s7(Number_type(12, 1), Number_type(10, 1));
  Point_2 t7(Number_type(14, 1), Number_type(10, 1));
  Arc_x_monotone_section_2 xc2(circ_2, s7, t7, CGAL::CLOCKWISE);
  x_curves.push_back(xc2);

  //create x-monotone polycurve
  X_monotone_polycurve x_polycurve_1 =
    traits.construct_x_monotone_curve_2_object()(x_curves.begin(),
                                                 x_curves.end());

  // Insert polycurves to Arangment and print.
  Polycurve_circ_arc_arrangment polycurve_arrangment(&traits);
  insert(polycurve_arrangment, polycurve_1);
  insert(polycurve_arrangment, polycurve_2);
  insert(polycurve_arrangment, x_polycurve_1);
  std::cout << "Arrangment Statistics: " << std::endl;
  print_arrangement(polycurve_arrangment);

  return 0;
}
#endif
