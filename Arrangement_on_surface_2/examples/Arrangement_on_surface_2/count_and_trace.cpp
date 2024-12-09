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

#include "arr_print.h"

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using Segment_base_traits = CGAL::Arr_segment_traits_2<Kernel>;
using Segment_cnt_traits = CGAL::Arr_counting_traits_2<Segment_base_traits>;
using Segment_traits = CGAL::Arr_tracing_traits_2<Segment_cnt_traits>;
using Segment_arrangement = CGAL::Arrangement_2<Segment_traits>;
using Point = Segment_traits::Point_2;
using Segment = Segment_traits::Curve_2;

using Geodesic_base_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
using Geodesic_cnt_traits = CGAL::Arr_counting_traits_2<Geodesic_base_traits>;
using Geodesic_traits = CGAL::Arr_tracing_traits_2<Geodesic_cnt_traits>;
using Topol_traits = CGAL::Arr_spherical_topology_traits_2<Geodesic_traits>;
using Geodesic_arrangement =
  CGAL::Arrangement_on_surface_2<Geodesic_traits, Topol_traits>;
using Geodesic_point = Geodesic_traits::Point_2;
using Geodesic_curve = Geodesic_traits::Curve_2;

using Nt_traits = CGAL::CORE_algebraic_number_traits;
using NT = Nt_traits::Rational;
using Rational = Nt_traits::Rational;
using Algebraic = Nt_traits::Algebraic;
using Rat_kernel = CGAL::Cartesian<Rational>;
using Alg_kernel = CGAL::Cartesian<Algebraic>;
using Rat_point = Rat_kernel::Point_2;
using Bezier_base_traits =
  CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
using Bezier_cnt_traits = CGAL::Arr_counting_traits_2<Bezier_base_traits>;
using Bezier_traits = CGAL::Arr_tracing_traits_2<Bezier_cnt_traits>;
using Bezier_arrangement = CGAL::Arrangement_2<Bezier_traits>;

int main() {
  Segment_traits seg_traits;
  seg_traits.disable_all_traces();
  seg_traits.enable_trace(Segment_traits::COMPARE_Y_AT_X_2_OP);
  Segment_arrangement seg_arr(&seg_traits);
  std::vector<Point> ps = { Point(0,0), Point(1,0), Point(0,1) };
  std::vector<Segment> segs =
    { Segment(ps[0], ps[1]), Segment(ps[1], ps[2]), Segment(ps[2], ps[0]) };
  CGAL::insert(seg_arr, segs.begin(), segs.end());
  std::cout << seg_traits;
  print_arrangement_size(seg_arr);
  std::cout << std::endl;

  Geodesic_traits geodesic_traits;
  geodesic_traits.disable_all_traces();
  geodesic_traits.enable_trace(Geodesic_traits::COMPARE_XY_2_OP);
  auto ctr_p = geodesic_traits.construct_point_2_object();
  auto ctr_cv = geodesic_traits.construct_curve_2_object();
  std::vector<Geodesic_point> gps =
    { ctr_p(-1,0,0), ctr_p(0,-1,0), ctr_p(0,0,-1) };
  std::vector<Geodesic_curve> gas =
    { ctr_cv(gps[0], gps[1]), ctr_cv(gps[1], gps[2]), ctr_cv(gps[2], gps[0]) };
  Geodesic_arrangement geodesic_arr(&geodesic_traits);
  CGAL::insert(geodesic_arr, gas.begin(), gas.end());
  std::cout << geodesic_traits;
  print_arrangement_size(geodesic_arr);

  Bezier_traits bezier_traits;
  bezier_traits.disable_all_traces();
  std::cout << bezier_traits;

  return 0;
}
