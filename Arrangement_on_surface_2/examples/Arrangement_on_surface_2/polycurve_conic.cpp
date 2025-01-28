// Testing the do_equal function

#include <CGAL/config.h>

#ifndef CGAL_USE_CORE

#include <iostream>

int main() {
  std::cout << "Sorry, this example needs CORE ...\n";
  return 0;
}

#else

#include <vector>

#include <CGAL/basic.h>
#include <CGAL/Arr_polycurve_traits_2.h>

#include "arr_conics.h"
#include "arr_print.h"

using Polycurve_conic_traits_2 = CGAL::Arr_polycurve_traits_2<Traits>;
using X_monotone_polycurve = Polycurve_conic_traits_2::X_monotone_curve_2;
using Polycurve = Polycurve_conic_traits_2::Curve_2;
using Polycurve_conic_arrangment = CGAL::Arrangement_2<Polycurve_conic_traits_2>;

int main() {
  Traits sub_traits;
  Polycurve_conic_traits_2 traits(&sub_traits);
  auto ctr_sub_cv = sub_traits.construct_curve_2_object();
  auto ctr_sub_xcv = sub_traits.construct_x_monotone_curve_2_object();

  // Polycurve construction functors
  auto ctr_xpolycurve = traits.construct_x_monotone_curve_2_object();
  auto ctr_polycurve = traits.construct_curve_2_object();

  // Containers to store conic curves that will be used to create polycurve.
  std::vector<Conic_arc> conic_curves;
  std::vector<X_monotone_conic_arc> xmono_conic_curves_2;

  // Create polycurves
  // y=x^2
  Conic_arc c3 = ctr_sub_cv(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                            Point(Algebraic(0), Algebraic(0)),
                            Point(Algebraic(3), Algebraic(9)));
  Conic_arc c4 = ctr_sub_cv(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                            Point(Algebraic(3), Algebraic(9)),
                            Point(Algebraic(5), Algebraic(25)));
  Conic_arc c5 = ctr_sub_cv(0,1,0,1,0,0, CGAL::COUNTERCLOCKWISE,
                            Point(Algebraic(-25), Algebraic(-5)),
                            Point(Algebraic(0), Algebraic(0)));

  Conic_arc c6 = ctr_sub_cv(1,1,0,6,-26,162,CGAL::COUNTERCLOCKWISE,
                            Point(Algebraic(-7), Algebraic(13)),
                            Point(Algebraic(-3), Algebraic(9)));
  Conic_arc c7 = ctr_sub_cv(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                            Point(Algebraic(-3), Algebraic(9)),
                            Point(Algebraic(0), Algebraic(0)));
  Conic_arc c8 = ctr_sub_cv(0,1,0,-1,0,0, CGAL::COUNTERCLOCKWISE,
                            Point(Algebraic(0), Algebraic(0)),
                            Point(Algebraic(4), Algebraic(-2)));

  Conic_arc c9 = ctr_sub_cv(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                            Point(Algebraic(-5), Algebraic(25)),
                            Point(Algebraic(5), Algebraic(25)));

  // Construct poly-curve
  conic_curves.clear();
  conic_curves.push_back(c9);
  Polycurve conic_polycurve_1 =
    ctr_polycurve(conic_curves.begin(), conic_curves.end());

  Conic_arc c11 = ctr_sub_cv(0,1,0,-1,0,0,CGAL::COUNTERCLOCKWISE,
                             Point(Algebraic(25), Algebraic(-5)),
                             Point(Algebraic(0), Algebraic(0)));
  Conic_arc c12 = ctr_sub_cv(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                             Point(Algebraic(0), Algebraic(0)),
                             Point(Algebraic(5), Algebraic(25)));

  // Construct poly-curve
  conic_curves.clear();
  conic_curves.push_back(c11);
  conic_curves.push_back(c12);
  Polycurve conic_polycurve_2 =
    ctr_polycurve(conic_curves.begin(), conic_curves.end());

  // Construct x-monotone conic curves from conic curves
  X_monotone_conic_arc xc3 = ctr_sub_xcv(c3);
  X_monotone_conic_arc xc4 = ctr_sub_xcv(c4);
  X_monotone_conic_arc xc5 = ctr_sub_xcv(c5);
  X_monotone_conic_arc xc6 = ctr_sub_xcv(c6);
  X_monotone_conic_arc xc7 = ctr_sub_xcv(c7);
  X_monotone_conic_arc xc8 = ctr_sub_xcv(c8);

  // Construct x-monotone poly-curve from x-monotone conic curves.
  xmono_conic_curves_2.clear();
  xmono_conic_curves_2.push_back(xc5);
  xmono_conic_curves_2.push_back(xc3);
  xmono_conic_curves_2.push_back(xc4);
  X_monotone_polycurve conic_x_mono_polycurve_1 =
    ctr_xpolycurve(xmono_conic_curves_2.begin(),
                   xmono_conic_curves_2.end());

  // Construct x-monotone poly-curve.
  xmono_conic_curves_2.clear();
  xmono_conic_curves_2.push_back(xc6);
  xmono_conic_curves_2.push_back(xc7);
  xmono_conic_curves_2.push_back(xc8);
  X_monotone_polycurve conic_x_mono_polycurve_2 =
    ctr_xpolycurve(xmono_conic_curves_2.begin(), xmono_conic_curves_2.end());

  // Insert the Polycurves into arrangement and print.
  Polycurve_conic_arrangment x_pc_arrangment(&traits);
  insert(x_pc_arrangment, conic_x_mono_polycurve_1);
  insert(x_pc_arrangment, conic_x_mono_polycurve_2);
  std::cout << "X-monotone polycurve arrangement Statistics:\n";
  print_arrangement(x_pc_arrangment);

  Polycurve_conic_arrangment pc_arrangment(&traits);
  insert(pc_arrangment, conic_polycurve_1);
  insert(pc_arrangment, conic_polycurve_2);
  std::cout << "Polycurve arrangement Statistics:\n";
  print_arrangement(pc_arrangment);

  return 0;
}

#endif
