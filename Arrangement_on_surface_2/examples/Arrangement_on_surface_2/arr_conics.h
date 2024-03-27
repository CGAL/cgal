#ifndef ARR_CONICS_H
#define ARR_CONICS_H

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>

using Nt_traits = CGAL::CORE_algebraic_number_traits;
using Rational = Nt_traits::Rational;
using Rat_kernel = CGAL::Cartesian<Rational>;
using Rat_point = Rat_kernel::Point_2;
using Rat_segment = Rat_kernel::Segment_2;
using Rat_circle = Rat_kernel::Circle_2;
using Algebraic = Nt_traits::Algebraic;
using Alg_kernel = CGAL::Cartesian<Algebraic>;
using Traits = CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
using Point = Traits::Point_2;
using Conic_arc = Traits::Curve_2;
using X_monotone_conic_arc = Traits::X_monotone_curve_2;
using Arrangement = CGAL::Arrangement_2<Traits>;

#endif
