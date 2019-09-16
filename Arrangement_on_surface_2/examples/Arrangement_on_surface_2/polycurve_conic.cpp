// Testing the do_equal function

#include <CGAL/config.h>

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
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include "arr_print.h"

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
  Conic_traits_2;
typedef Conic_traits_2::Point_2                       Conic_point_2;
typedef Conic_traits_2::Curve_2                       Conic_curve_2;
typedef Conic_traits_2::X_monotone_curve_2            Conic_x_monotone_curve_2;
typedef CGAL::Arr_polycurve_traits_2<Conic_traits_2>  Polycurve_conic_traits_2;
typedef Polycurve_conic_traits_2::X_monotone_curve_2  X_monotone_polycurve;
typedef Polycurve_conic_traits_2::Curve_2             Polycurve;
typedef CGAL::Arrangement_2<Polycurve_conic_traits_2> Polycurve_conic_arrangment;

int main()
{
  Polycurve_conic_traits_2 traits;

  // Polycurve construction functors
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2
    construct_x_mono_polycurve = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Construct_curve_2  construct_polycurve =
    traits.construct_curve_2_object();

  // Containers to store conic curves that will be used to create polycurve.
  std::vector<Conic_curve_2> conic_curves;
  std::vector<Conic_x_monotone_curve_2> xmono_conic_curves_2;

  // Create polycurves
  // y=x^2
  Conic_curve_2 c3(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(0), Algebraic(0)),
                   Conic_point_2(Algebraic(3), Algebraic(9)));
  Conic_curve_2 c4(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(3), Algebraic(9)),
                   Conic_point_2(Algebraic(5), Algebraic(25)));
  Conic_curve_2 c5(0,1,0,1,0,0, CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(-25), Algebraic(-5)),
                   Conic_point_2(Algebraic(0), Algebraic(0)));

  Conic_curve_2 c6(1,1,0,6,-26,162,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(-7), Algebraic(13)),
                   Conic_point_2(Algebraic(-3), Algebraic(9)));
  Conic_curve_2 c7(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(-3), Algebraic(9)),
                   Conic_point_2(Algebraic(0), Algebraic(0)));
  Conic_curve_2 c8(0,1,0,-1,0,0, CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(0), Algebraic(0)),
                   Conic_point_2(Algebraic(4), Algebraic(-2)));

  Conic_curve_2 c9(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(-5), Algebraic(25)),
                   Conic_point_2(Algebraic(5), Algebraic(25)));

  // Construct poly-curve
  conic_curves.clear();
  conic_curves.push_back(c9);
  Polycurve conic_polycurve_1 =
    construct_polycurve(conic_curves.begin(), conic_curves.end());

  Conic_curve_2 c11(0,1,0,-1,0,0,CGAL::COUNTERCLOCKWISE,
                    Conic_point_2(Algebraic(25), Algebraic(-5)),
                    Conic_point_2(Algebraic(0), Algebraic(0)));
  Conic_curve_2 c12(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                    Conic_point_2(Algebraic(0), Algebraic(0)),
                    Conic_point_2(Algebraic(5), Algebraic(25)));

  // Construct poly-curve
  conic_curves.clear();
  conic_curves.push_back(c11);
  conic_curves.push_back(c12);
  Polycurve conic_polycurve_2 =
    construct_polycurve(conic_curves.begin(), conic_curves.end());

  // Construct x-monotone conic curves from conic curves
  Conic_x_monotone_curve_2 xc3(c3);
  Conic_x_monotone_curve_2 xc4(c4);
  Conic_x_monotone_curve_2 xc5(c5);
  Conic_x_monotone_curve_2 xc6(c6);
  Conic_x_monotone_curve_2 xc7(c7);
  Conic_x_monotone_curve_2 xc8(c8);

  // Construct x-monotone poly-curve from x-monotone conic curves.
  xmono_conic_curves_2.clear();
  xmono_conic_curves_2.push_back(xc5);
  xmono_conic_curves_2.push_back(xc3);
  xmono_conic_curves_2.push_back(xc4);
  X_monotone_polycurve conic_x_mono_polycurve_1 =
    construct_x_mono_polycurve(xmono_conic_curves_2.begin(),
                               xmono_conic_curves_2.end());

  // Construct x-monotone poly-curve.
  xmono_conic_curves_2.clear();
  xmono_conic_curves_2.push_back(xc6);
  xmono_conic_curves_2.push_back(xc7);
  xmono_conic_curves_2.push_back(xc8);
  X_monotone_polycurve conic_x_mono_polycurve_2 =
    construct_x_mono_polycurve(xmono_conic_curves_2.begin(),
                               xmono_conic_curves_2.end());

  // Insert the Polycurves into arrangment and print.
  Polycurve_conic_arrangment x_pc_arrangment(&traits);
  insert(x_pc_arrangment, conic_x_mono_polycurve_1);
  insert(x_pc_arrangment, conic_x_mono_polycurve_2);
  std::cout << "X-monotone polycurve arrangement Statistics: " << std::endl;
  print_arrangement(x_pc_arrangment);

  Polycurve_conic_arrangment pc_arrangment(&traits);
  insert(pc_arrangment, conic_polycurve_1);
  insert(pc_arrangment, conic_polycurve_2);
  std::cout << "Polycurve arrangement Statistics: " << std::endl;
  print_arrangement(pc_arrangment);

  return 0;
}

#endif
