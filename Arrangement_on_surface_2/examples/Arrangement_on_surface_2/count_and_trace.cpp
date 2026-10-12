#include <iostream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_counting_traits_2.h>
#include <CGAL/Arr_tracing_traits_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Arr_rational_function_traits_2.h>

#include "arr_print.h"

// #define COUNT_TRACE_ORDER 1

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

// Segments
using Segment_base_traits = CGAL::Arr_segment_traits_2<Kernel>;
#if COUNT_TRACE_ORDER
using Segment_cnt_traits = CGAL::Arr_counting_traits_2<Segment_base_traits>;
using Segment_trc_traits = CGAL::Arr_tracing_traits_2<Segment_cnt_traits>;
using Segment_traits = Segment_trc_traits;
#else
using Segment_trc_traits = CGAL::Arr_tracing_traits_2<Segment_base_traits>;
using Segment_cnt_traits = CGAL::Arr_counting_traits_2<Segment_trc_traits>;
using Segment_traits = Segment_cnt_traits;
#endif
using Segment_arrangement = CGAL::Arrangement_2<Segment_traits>;
using Point = Segment_traits::Point_2;
using Segment = Segment_traits::Curve_2;

// Geodesic arcs
using Geodesic_base_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
#if COUNT_TRACE_ORDER
using Geodesic_cnt_traits = CGAL::Arr_counting_traits_2<Geodesic_base_traits>;
using Geodesic_trc_traits = CGAL::Arr_tracing_traits_2<Geodesic_cnt_traits>;
using Geodesic_traits = Geodesic_trc_traits;
#else
using Geodesic_trc_traits = CGAL::Arr_tracing_traits_2<Geodesic_base_traits>;
using Geodesic_cnt_traits = CGAL::Arr_counting_traits_2<Geodesic_trc_traits>;
using Geodesic_traits = Geodesic_cnt_traits;
#endif
using Topol_traits = CGAL::Arr_spherical_topology_traits_2<Geodesic_traits>;
using Geodesic_arrangement = CGAL::Arrangement_on_surface_2<Geodesic_traits, Topol_traits>;
using Geodesic_point = Geodesic_traits::Point_2;
using Geodesic_curve = Geodesic_traits::Curve_2;

// Bezier curves
using Nt_traits = CGAL::CORE_algebraic_number_traits;
using NT = Nt_traits::Rational;
using Rational = Nt_traits::Rational;
using Algebraic = Nt_traits::Algebraic;
using Rat_kernel = CGAL::Cartesian<Rational>;
using Alg_kernel = CGAL::Cartesian<Algebraic>;
using Rat_point = Rat_kernel::Point_2;
using Bezier_base_traits = CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
#if COUNT_TRACE_ORDER
using Bezier_cnt_traits = CGAL::Arr_counting_traits_2<Bezier_base_traits>;
using Bezier_trc_traits = CGAL::Arr_tracing_traits_2<Bezier_cnt_traits>;
using Bezier_traits = Bezier_trc_traits;
#else
using Bezier_trc_traits = CGAL::Arr_tracing_traits_2<Bezier_base_traits>;
using Bezier_cnt_traits = CGAL::Arr_counting_traits_2<Bezier_trc_traits>;
using Bezier_traits = Bezier_cnt_traits;
#endif
using Bezier_arrangement = CGAL::Arrangement_2<Bezier_traits>;

// Rational function arcs
using Rat_func_number_type = CORE::BigInt;
using Rat_func_algebraic_kernel = CGAL::Algebraic_kernel_d_1<Rat_func_number_type>;
using Rat_func_base_traits = CGAL::Arr_rational_function_traits_2<Rat_func_algebraic_kernel>;
using Rat_func_cnt_traits = CGAL::Arr_counting_traits_2<Rat_func_base_traits>;
using Rat_func_trc_traits = CGAL::Arr_tracing_traits_2<Rat_func_base_traits>;
#if COUNT_TRACE_ORDER
using Rat_func_traits = CGAL::Arr_tracing_traits_2<Bezier_cnt_traits>;
#else
using Rat_func_traits = CGAL::Arr_counting_traits_2<Rat_func_trc_traits>;
#endif
using Rat_func_arrangement = CGAL::Arrangement_2<Rat_func_traits>;

int main() {
  // Count and trace segment traits
  Segment_traits seg_traits;
#if COUNT_TRACE_ORDER
  seg_traits.disable_all_traces();
  seg_traits.enable_trace(Segment_trc_traits::COMPARE_Y_AT_X_2_OP);
#else
  seg_traits.traits().disable_all_traces();
  seg_traits.traits().enable_trace(Segment_trc_traits::COMPARE_Y_AT_X_2_OP);
#endif
  Segment_arrangement seg_arr(&seg_traits);
  std::vector<Point> ps = { Point(0,0), Point(1,0), Point(0,1) };
  std::vector<Segment> segs = { Segment(ps[0], ps[1]), Segment(ps[1], ps[2]), Segment(ps[2], ps[0]) };
  CGAL::insert(seg_arr, segs.begin(), segs.end());
#if COUNT_TRACE_ORDER
  std::cout << seg_traits.traits();
#else
  std::cout << seg_traits;
#endif
  print_arrangement_size(seg_arr);
  std::cout << std::endl;

  // Count and trace geodesic arc traits
  Geodesic_traits geodesic_traits;
#if COUNT_TRACE_ORDER
  geodesic_traits.disable_all_traces();
  geodesic_traits.enable_trace(Geodesic_trc_traits::COMPARE_XY_2_OP);
#else
  geodesic_traits.traits().disable_all_traces();
  geodesic_traits.traits().enable_trace(Geodesic_trc_traits::COMPARE_XY_2_OP);
#endif
  auto ctr_p = geodesic_traits.construct_point_2_object();
  auto ctr_cv = geodesic_traits.construct_curve_2_object();
  std::vector<Geodesic_point> gps = { ctr_p(-1,0,0), ctr_p(0,-1,0), ctr_p(0,0,-1) };
  std::vector<Geodesic_curve> gas = { ctr_cv(gps[0], gps[1]), ctr_cv(gps[1], gps[2]), ctr_cv(gps[2], gps[0]) };
  Geodesic_arrangement geodesic_arr(&geodesic_traits);
  CGAL::insert(geodesic_arr, gas.begin(), gas.end());
#if COUNT_TRACE_ORDER
  std::cout << geodesic_traits.traits();
#else
  std::cout << geodesic_traits;
#endif
  print_arrangement_size(geodesic_arr);

  // Count and trace Bezier curve traits
  Bezier_traits bezier_traits;
#if COUNT_TRACE_ORDER
  bezier_traits.disable_all_traces();
  std::cout << bezier_traits.traits();
#else
  bezier_traits.traits().disable_all_traces();
  std::cout << bezier_traits;
#endif

  // Count and trace rational function traits
  auto rat_func_base_traits = std::shared_ptr<Rat_func_base_traits>(new Rat_func_base_traits);
#if COUNT_TRACE_ORDER
  auto rat_func_mid_traits = std::shared_ptr<Rat_func_cnt_traits>(new Rat_func_cnt_traits(rat_func_base_traits));
#else
  auto rat_func_mid_traits = std::shared_ptr<Rat_func_trc_traits>(new Rat_func_trc_traits(rat_func_base_traits));
#endif
  auto rat_func_traits = std::shared_ptr<Rat_func_traits>(new Rat_func_traits(rat_func_mid_traits));
#if COUNT_TRACE_ORDER
  rat_func_traits->disable_all_traces();
  rat_func_traits->enable_trace(Rat_func_trc_traits::CONSTRUCT_POINT_2_OP);
  rat_func_traits->enable_trace(Rat_func_trc_traits::CONSTRUCT_POINT_2_XY_OP);
#else
  rat_func_traits->shared_traits()->disable_all_traces();
  rat_func_traits->shared_traits()->enable_trace(Rat_func_trc_traits::CONSTRUCT_POINT_2_OP);
  rat_func_traits->shared_traits()->enable_trace(Rat_func_trc_traits::CONSTRUCT_POINT_2_XY_OP);
#endif
  rat_func_traits->construct_point_2_object()(0, 0);
#if COUNT_TRACE_ORDER
  std::cout << *(rat_func_traits->shared_traits();
#else
  std::cout << *rat_func_traits;
#endif

  return 0;
}
