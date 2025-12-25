#ifndef CGAL_SURFACE_SWEEP_2_DEFAULT_ARR_TRAITS_H
#define CGAL_SURFACE_SWEEP_2_DEFAULT_ARR_TRAITS_H

#include <CGAL/license/Surface_sweep_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_rational_function_traits_2.h>

/*! \file
 *
 * Definaition of the Default_arr_traits.
 */

namespace CGAL {
namespace Surface_sweep_2 {

//!
template <typename Curve> struct Default_arr_traits {};

//!
template <typename Kernel>
struct Default_arr_traits<CGAL::Segment_2<Kernel>> {
  using Traits = CGAL::Arr_segment_traits_2<Kernel>;
};

//!
template <typename Kernel>
struct Default_arr_traits<CGAL::Arr_segment_2<Kernel>> {
  using Traits = CGAL::Arr_segment_traits_2<Kernel>;
};

//!
template <typename Kernel>
struct Default_arr_traits<CGAL::internal::Polycurve_2<CGAL::Arr_segment_2<Kernel>, typename Kernel::Point_2>> {
  using Subtraits = CGAL::Arr_segment_traits_2<Kernel>;
  using Traits = CGAL::Arr_polyline_traits_2<Subtraits>;
};

//!
template <typename Rat_kernel_, class Alg_kernel_, class Nt_traits_>
struct Default_arr_traits<CGAL::Conic_arc_2<Rat_kernel_, Alg_kernel_, Nt_traits_>> {
  using Traits = CGAL::Arr_conic_traits_2<Rat_kernel_, Alg_kernel_, Nt_traits_>;
};

//!
// template <typename AlgebraicKernel_d_1>
// class Arr_rational_function_traits_2;

// namespace Arr_rational_arc {

// //!
// template <typename Algebraic_kernel_>
// class Rational_arc_d_1;

// }

//!
template <typename Algebraic_kernel_>
struct Default_arr_traits<CGAL::Arr_rational_arc::Rational_arc_d_1<Algebraic_kernel_>> {
  using Traits = CGAL::Arr_rational_function_traits_2<Algebraic_kernel_>;
};

//!
template <typename Kernel_, bool Filter_>
struct Default_arr_traits<CGAL::_Circle_segment_2<Kernel_, Filter_>> {
  using Traits = CGAL::Arr_circle_segment_traits_2<Kernel_, Filter_>;
};

//!
template <typename Kernel>
struct Default_arr_traits<CGAL::Arr_linear_object_2<Kernel>> {
  using Traits = CGAL::Arr_linear_traits_2<Kernel>;
};

}
}

#endif
