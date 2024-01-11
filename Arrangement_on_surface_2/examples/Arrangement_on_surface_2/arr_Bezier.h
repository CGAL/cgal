#ifndef ARR_BEZIER_H
#define ARR_BEZIER_H

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arrangement_2.h>

using Nt_traits = CGAL::CORE_algebraic_number_traits;
using NT = Nt_traits::Rational;
using Rational = Nt_traits::Rational;
using Algebraic = Nt_traits::Algebraic;
using Rat_kernel = CGAL::Cartesian<Rational>;
using Alg_kernel = CGAL::Cartesian<Algebraic>;
using Rat_point = Rat_kernel::Point_2;
using Traits =
  CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
using Bezier_x_monotone_curve = Traits::X_monotone_curve_2;
using Bezier_curve = Traits::Curve_2;
using Arrangement = CGAL::Arrangement_2<Traits>;

#endif
