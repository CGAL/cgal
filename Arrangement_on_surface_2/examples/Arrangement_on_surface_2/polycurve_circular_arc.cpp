
// Constructing an arrangement of polycurves.

#include <CGAL/basic.h>
#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
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
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include "arr_print.h"

///////////////
//circle segment traits
//////////////
typedef CGAL::Quotient<CGAL::MP_Float>                                    Number_type;
typedef CGAL::Cartesian<Number_type>                                      Kernel;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>                         Arc_traits_2;
typedef CGAL::Arr_polyline_traits_2<Arc_traits_2>                         Polycurve_arc_traits_2;

typedef Arc_traits_2::CoordNT                                             CoordNT;
typedef Arc_traits_2::Point_2                                             Point_2;
typedef Arc_traits_2::Curve_2                                             Arc_section_2;
typedef Arc_traits_2::X_monotone_curve_2                                  Arc_x_monotone_section_2;

typedef Polycurve_arc_traits_2::X_monotone_curve_2                        X_monotone_polycurve;
typedef Polycurve_arc_traits_2::Curve_2                                   Polycurve;
typedef Kernel::Circle_2                                                  Circle_2;
typedef CGAL::Arrangement_2<Polycurve_arc_traits_2>                       Polycurve_circ_arc_arrangment;


int main ()
{

  Polycurve_arc_traits_2 traits;

  //Containers to store conic curves that will be used to create polycurve.
  std::vector<Arc_section_2> curves;
  std::vector<Arc_x_monotone_section_2> x_curves;

  // Create a circular arc of the circle, directed clockwise from
  // (-1, 0) to (1, 0). Note that we orient the
  // supporting circle accordingly.
  Kernel::Point_2  c6 = Kernel::Point_2 (0, 0);
  Point_2          s6 = Point_2 ( Number_type (-1, 1), Number_type (0, 1) );
  Point_2          t6 = Point_2 ( Number_type (1, 1), Number_type (0, 1));
  Arc_section_2    circ_arc1 (c6, 1, CGAL::CLOCKWISE, s6, t6);
  curves.push_back (circ_arc1);

  // Create a circular arc of the unit circle, directed clockwise from
  // (1, 0) to (3, 0). Note that we orient the
  // supporting circle accordingly.
  Kernel::Point_2  c1 = Kernel::Point_2 (3, 0);
  Point_2          s1 = Point_2 ( Number_type (1, 1), Number_type (0, 1) );
  Point_2          t1 = Point_2 ( Number_type (5, 1), Number_type (0, 1));
  Arc_section_2    circ_arc2 (c1, 2, CGAL::CLOCKWISE, s1, t1);
  curves.push_back (circ_arc2);

  //create polycurve
  Polycurve polycurve_1 = traits.construct_curve_2_object()(curves.begin(), curves.end() );

  //Another polycurve
  curves.clear();
  
  Kernel::Point_2  center = Kernel::Point_2 (-10, 10);
  Point_2          source = Point_2 ( Number_type (-10, 1), Number_type (13, 1) );
  Point_2          target = Point_2 ( Number_type (-7, 1), Number_type (10, 1));
  Arc_section_2    circ_arc (center, 3, CGAL::CLOCKWISE, source, target);
  curves.push_back( circ_arc );

  center = Kernel::Point_2 (-20, 10);
  source = Point_2 ( Number_type (-7, 1), Number_type (10, 1) );
  target = Point_2 ( Number_type (-20, 1), Number_type (23, 1));
  circ_arc =  Arc_section_2(center, 13, CGAL::CLOCKWISE, source, target);
  curves.push_back( circ_arc );

  center = Kernel::Point_2 (-20, 25);
  source = Point_2 ( Number_type (-20, 1), Number_type (23, 1) );
  target = Point_2 ( Number_type (-20, 1), Number_type (27, 1));
  circ_arc =  Arc_section_2(center, 2, CGAL::CLOCKWISE, source, target);
  curves.push_back( circ_arc );

  Polycurve polycurve_2 = traits.construct_curve_2_object()(curves.begin(), curves.end() );

  //circle to be used by x-monotone polycurve
  Circle_2 circ = Circle_2(c6, 1, CGAL::CLOCKWISE);
  Arc_x_monotone_section_2 xc1(circ, s6, t6, CGAL::CLOCKWISE);
  x_curves.push_back(xc1);

  Circle_2 circ1 = Circle_2(c1, 2, CGAL::CLOCKWISE);
  Arc_x_monotone_section_2 xc2(circ1, s1, t1, CGAL::CLOCKWISE);
  x_curves.push_back(xc2);

  //create x-monotone polycurve
  X_monotone_polycurve x_polycurve_1 = traits.construct_x_monotone_curve_2_object()( x_curves.begin(), x_curves.end() );

  Kernel::Point_2  cen = Kernel::Point_2 (-2, 0);
  Circle_2 circ2 = Circle_2(cen, 2, CGAL::CLOCKWISE);
  Point_2          s2 = Point_2 ( Number_type (-4, 1), Number_type (0, 1) );
  Point_2          t2 = Point_2 ( Number_type (0, 1), Number_type (0, 1));
  Arc_x_monotone_section_2 xc3(circ2, s2, t2, CGAL::CLOCKWISE);
  
  Arc_x_monotone_section_2 xc4( Kernel::Point_2(0, 0), Kernel::Point_2(2, 0) );

  //create x-monotone polycurve
  x_curves.clear();

  
  x_curves.push_back(xc3);
  x_curves.push_back(xc4);

  X_monotone_polycurve x_polycurve_2 = traits.construct_x_monotone_curve_2_object()( x_curves.begin(), x_curves.end() );

  //Print the Polycurves.
  std::cout << std::endl;
  std::cout << "Polycurve 1: " << polycurve_1 << std::endl;
  std::cout << "Polycurve 2: " << polycurve_2 << std::endl;
  std::cout << "X-monotone Polycurve 1: " << x_polycurve_1 << std::endl;
  std::cout << "X-monotone Polycurve 2: " << x_polycurve_2 << std::endl;

  //create opposite of conic_x_mono_polycurve_2
  X_monotone_polycurve opposite_x_monotone_polycurve_1 = traits.construct_opposite_2_object()(x_polycurve_1);
  //std::cout << "Opposite of X-monotone Polycurve 1: " << opposite_x_monotone_polycurve_1 << std::endl;

  //Waqar: fix this once the intersect functor for circular arc polycurve is fixed. 
  Polycurve_circ_arc_arrangment polycurve_arrangment(&traits);
  insert(polycurve_arrangment, polycurve_1);
  insert(polycurve_arrangment, polycurve_2);
  insert(polycurve_arrangment, x_polycurve_1);
  insert(polycurve_arrangment, x_polycurve_2);
  std::cout << "Arrangment Statistics: " << std::endl; 
   print_arrangement (polycurve_arrangment);


  return 0;
}
#endif
