
#include "CGAL/Arr_algebraic_segment_traits_2.h"
#include "CGAL/Arr_circle_segment_traits_2.h"
#include "CGAL/Arr_conic_traits_2.h"
#include "CGAL/Arr_default_dcel.h"
#include "CGAL/Arr_enums.h"
#include "CGAL/Arr_linear_traits_2.h"
#include "CGAL/Arr_rational_function_traits_2.h"
#include "CGAL/Arr_segment_traits_2.h"
#include "CGAL/Arr_spherical_topology_traits_2.h"
#include "CGAL/Arrangement_2.h"
#include "CGAL/Arrangement_on_surface_2.h"
#include "CGAL/CORE_algebraic_number_traits.h"
#include "CGAL/Draw_aos/Arr_viewer.h"
#include <array>
#include <fstream>
#include <iostream>
#include "CGAL/Exact_predicates_exact_constructions_kernel.h"
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include <CGAL/draw_arrangement_2.h>

// void draw_segments_arr_1() {
//   using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
//   using Segment_traits = CGAL::Arr_segment_traits_2<Exact_kernel>;
//   using Point_2 = Segment_traits::Point_2;
//   using Arrangement = CGAL::Arrangement_2<Segment_traits>;
//   // Make a square
//   Arrangement arr;
//   auto traits = arr.traits();
//   auto cst_x_curve = traits->construct_x_monotone_curve_2_object();
//   auto square = {cst_x_curve({0, 0}, {5, 0}), cst_x_curve({5, 0}, {5, 5}), cst_x_curve({5, 5}, {0, 5}),
//                  cst_x_curve({0, 5}, {0, 0})};
//   insert(arr, square.begin(), square.end());
//   auto hole_triangle = {
//       cst_x_curve({1, 1}, {2, 1}),
//       cst_x_curve({2, 1}, {2, 2}),
//       cst_x_curve({2, 2}, {1, 1}),
//   };
//   insert(arr, hole_triangle.begin(), hole_triangle.end());
//   // CGAL::draw_viewer(arr);
// }

// void draw_segments_arr_2() {
//   using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
//   using Segment_traits = CGAL::Arr_segment_traits_2<Exact_kernel>;
//   using Point_2 = Segment_traits::Point_2;
//   using X_monotone_curve_2 = Segment_traits::X_monotone_curve_2;
//   using Arrangement = CGAL::Arrangement_2<Segment_traits>;
//   Arrangement arr;
//   auto traits = arr.traits();
//   auto cst_x_curve = traits->construct_x_monotone_curve_2_object();
//   // make a hexagon centered at the origin with radius 10
//   double radius = 10.0;
//   std::array<Point_2, 6> hexagon_points = {{
//       {radius * cos(0 * CGAL_PI / 3), radius * sin(0 * CGAL_PI / 3)}, // 0°
//       {radius * cos(1 * CGAL_PI / 3), radius * sin(1 * CGAL_PI / 3)}, // 60°
//       {radius * cos(2 * CGAL_PI / 3), radius * sin(2 * CGAL_PI / 3)}, // 120°
//       {radius * cos(3 * CGAL_PI / 3), radius * sin(3 * CGAL_PI / 3)}, // 180°
//       {radius * cos(4 * CGAL_PI / 3), radius * sin(4 * CGAL_PI / 3)}, // 240°
//       {radius * cos(5 * CGAL_PI / 3), radius * sin(5 * CGAL_PI / 3)}, // 300°
//   }};
//   std::array<X_monotone_curve_2, 6> hexagon;
//   for(size_t i = 0; i < hexagon_points.size(); ++i) {
//     size_t next_i = (i + 1) % hexagon_points.size();
//     hexagon[i] = cst_x_curve(hexagon_points[i], hexagon_points[next_i]);
//   }
//   // rect hole
//   auto hole_rectangle = {cst_x_curve({-2, -2}, {2, -2}), cst_x_curve({2, -2}, {2, 2}), cst_x_curve({2, 2}, {-2, 2}),
//                          cst_x_curve({-2, 2}, {-2, -2})};
//   // iso vertex inside rect hole
//   auto iso_vertex_inside_hole = Point_2{0.5, 0.5};
//   // degenerate segment below the rect hole
//   auto degenerate_segment = cst_x_curve({0, -3}, {1, -3});

//   CGAL::insert_point(arr, iso_vertex_inside_hole);
//   CGAL::insert(arr, hexagon.begin(), hexagon.end());
//   CGAL::insert(arr, hole_rectangle.begin(), hole_rectangle.end());
//   CGAL::insert(arr, degenerate_segment);
//   CGAL::draw_viewer(arr);
// }

// void draw_segments_arr_3() {
//   // generate random segments and draw them
//   using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
//   using Segment_traits = CGAL::Arr_segment_traits_2<Exact_kernel>;
//   using Point_2 = Segment_traits::Point_2;
//   using Arrangement = CGAL::Arrangement_2<Segment_traits>;
//   using X_monotone_curve_2 = Segment_traits::X_monotone_curve_2;
//   using Random = CGAL::Random;

//   Arrangement arr;
//   auto traits = arr.traits();
//   auto cst_x_curve = traits->construct_x_monotone_curve_2_object();
//   Random random;
//   std::vector<X_monotone_curve_2> segments;
//   std::cout << "Generating random segments..." << std::endl;
//   for(int i = 0; i < 100; ++i) {
//     // generate random points
//     Point_2 p1(random.get_double(-100, 100), random.get_double(-100, 100));
//     Point_2 p2(random.get_double(-100, 100), random.get_double(-100, 100));
//     // create a segment
//     X_monotone_curve_2 seg = cst_x_curve(p1, p2);
//     segments.push_back(seg);
//   }

//   std::cout << "Inserting segments into the arrangement..." << std::endl;
//   // insert segments into the arrangement
//   CGAL::insert(arr, segments.begin(), segments.end());

//   // draw the arrangement
//   // CGAL::draw_viewer(arr);
// }

void draw_linear_arr_1() {
  using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Traits = CGAL::Arr_linear_traits_2<Exact_kernel>;
  using Point_2 = Traits::Point_2;
  using Line_2 = Traits::Line_2;
  using Ray_2 = Traits::Ray_2;
  using Curve_2 = Traits::Curve_2;
  using Arrangement = CGAL::Arrangement_2<Traits>;
  using Face_const_handle = Arrangement::Face_const_handle;
  using Halfedge_const_handle = Arrangement::Halfedge_const_iterator;
  using X_monotone_curve_2 = Traits::X_monotone_curve_2;

  Arrangement arr;
  auto x_axis = X_monotone_curve_2(Ray_2(Point_2(0, 0), Point_2(1, 0)));
  auto y_axis = X_monotone_curve_2(Ray_2(Point_2(0, 0), Point_2(0, 1)));
  CGAL::insert(arr, x_axis);
  CGAL::insert(arr, y_axis);

  CGAL::draw_viewer(arr);
}

void draw_linear_arr_2() {
  using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Traits = CGAL::Arr_linear_traits_2<Exact_kernel>;
  using Point_2 = Traits::Point_2;
  using Line_2 = Traits::Line_2;
  using Segment_2 = Traits::Segment_2;
  using Ray_2 = Traits::Ray_2;
  using Curve_2 = Traits::Curve_2;
  using Arrangement = CGAL::Arrangement_2<Traits>;
  using Face_const_handle = Arrangement::Face_const_handle;
  using Halfedge_const_handle = Arrangement::Halfedge_const_iterator;
  using X_monotone_curve_2 = Traits::X_monotone_curve_2;

  Arrangement arr;
  auto& traits = *arr.traits();
  // Insert a n*n grid, each cell is a square of size 5
  int n = 1;
  // for(int i = 0; i < n; ++i) {
  //   Point_2 p1(i * 5, 0);
  //   Point_2 p2(i * 5, 1);
  //   CGAL::insert(arr, Curve_2(Line_2(p1, p2)));
  // }
  // for(int i = 0; i < n; ++i) {
  //   Point_2 p1(0, i * 5);
  //   Point_2 p2(1, i * 5);
  //   CGAL::insert(arr, Curve_2(Line_2(p1, p2)));
  // }
  // Generate a inner square(2*2) for all cells
  // And an inner triangle for each square
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      Point_2 p1(i * 5 + 1, j * 5 + 1);
      Point_2 p2(i * 5 + 4, j * 5 + 4);
      CGAL::insert(arr, Curve_2(Segment_2(p1, Point_2(p2.x(), p1.y()))));
      CGAL::insert(arr, Curve_2(Segment_2(Point_2(p1.x(), p2.y()), p2)));
      CGAL::insert(arr, Curve_2(Segment_2(p1, Point_2(p1.x(), p2.y()))));
      CGAL::insert(arr, Curve_2(Segment_2(Point_2(p2.x(), p1.y()), p2)));

      // Insert a triangle inside the square
      Point_2 tri_p1(i * 5 + 2, j * 5 + 2);
      Point_2 tri_p2(i * 5 + 3, j * 5 + 2);
      Point_2 tri_p3(i * 5 + 2.5, j * 5 + 3);
      CGAL::insert(arr, Curve_2(Segment_2(tri_p1, tri_p2)));
      CGAL::insert(arr, Curve_2(Segment_2(tri_p2, tri_p3)));
      CGAL::insert(arr, Curve_2(Segment_2(tri_p3, tri_p1)));

      // Connect the triangle to the square
      Point_2 top(i * 5 + 2.5, j * 5 + 4);
      CGAL::insert(arr, Curve_2(Segment_2(tri_p1, top)));
    }
  }

  auto arr2 = arr;

  CGAL::draw_viewer(arr);
}

// // supports segments
// void draw_circle_segs_arr() {
//   using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
//   using Traits = CGAL::Arr_circle_segment_traits_2<Exact_kernel>;
//   using Point_2 = Traits::Point_2;
//   using Curve_2 = Traits::Curve_2;
//   using Arrangement = CGAL::Arrangement_2<Traits>;

//   auto traits = Traits();
//   Arrangement arr;
//   auto cv1 = Curve_2(Exact_kernel::Circle_2({0, 0}, 10));
//   CGAL::insert(arr, cv1);
//   CGAL::draw(arr);
// }

// void draw_conic_arcs_arr() {
//   using Nt_traits = CGAL::CORE_algebraic_number_traits;
//   using Rational = Nt_traits::Rational;
//   using Rat_kernel = CGAL::Cartesian<Rational>;
//   using Rat_point = Rat_kernel::Point_2;
//   using Rat_segment = Rat_kernel::Segment_2;
//   using Rat_circle = Rat_kernel::Circle_2;
//   using Algebraic = Nt_traits::Algebraic;
//   using Alg_kernel = CGAL::Cartesian<Algebraic>;
//   using Traits = CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
//   using Point = Traits::Point_2;
//   using Conic_arc = Traits::Curve_2;
//   using X_monotone_conic_arc = Traits::X_monotone_curve_2;
//   using Arrangement = CGAL::Arrangement_2<Traits>;

//   Arrangement arr;
//   auto traits = Traits();
//   auto cst_x_curve = traits.construct_curve_2_object();

//   auto vert_seg = cst_x_curve(Rat_segment(Rat_point(0, 0), Rat_point(1, 0)));
//   auto hor_seg = cst_x_curve(Rat_segment(Rat_point(0, 0), Rat_point(0, 1)));

//   CGAL::insert(arr, vert_seg);
//   CGAL::insert(arr, hor_seg);

//   CGAL::draw(arr);
// }

// void draw_algebraic_arr() {
// #if CGAL_USE_GMP && CGAL_USE_MPFI
// #include <CGAL/Gmpz.h>
//   using Integer = CGAL::Gmpz;
// #elif CGAL_USE_CORE
// #include <CGAL/CORE_BigInt.h>
//   using Integer = CORE::BigInt;
// #else
// #include <CGAL/leda_integer.h>
//   using Integer = LEDA::integer;
// #endif
//   using Traits = CGAL::Arr_algebraic_segment_traits_2<Integer>;
//   using Arrangement = CGAL::Arrangement_2<Traits>;
//   using Polynomial = Traits::Polynomial_2;
//   using X_monotone_curve_2 = Traits::X_monotone_curve_2;
//   using Parameter_space_in_x_2 = Traits::Parameter_space_in_x_2;

//   Arrangement arr;
//   auto traits = arr.traits();
//   X_monotone_curve_2 cv;
//   auto param_space_in_x = traits->parameter_space_in_x_2_object();
//   auto ctr_cv = traits->construct_curve_2_object();
//   Polynomial x = CGAL::shift(Polynomial(1), 1, 0);
//   Polynomial y = CGAL::shift(Polynomial(1), 1, 1);
//   auto cst_x_curve = traits->construct_x_monotone_segment_2_object();
//   auto curve = ctr_cv(CGAL::ipower(x, 4) + CGAL::ipower(y, 3) - 1);
//   CGAL::insert(arr, curve);
//   // CGAL::draw(arr);
// }

// void draw_rational_arr() {
//   using AK1 = CGAL::Algebraic_kernel_d_1<CORE::BigInt>;
//   using Traits = CGAL::Arr_rational_function_traits_2<AK1>;
//   using Arrangement = CGAL::Arrangement_2<Traits>;
//   using Polynomial = Traits::Polynomial_1;
//   using Alg_real = Traits::Algebraic_real_1;
//   using Bound = Traits::Bound;

//   auto traits = Traits();
//   auto approx = traits.approximate_2_object();
//   auto cst_x_curve = traits.construct_x_monotone_curve_2_object();
//   Arrangement arr;
//   Polynomial x = CGAL::shift(Polynomial(1), 1);
//   Polynomial P1 = CGAL::ipower(x, 4) - 6 * x * x + 8;
//   Alg_real l(Bound(-2.1)), r(Bound(2.1));
//   auto cv1 = cst_x_curve(P1, l, r);
//   CGAL::insert(arr, cv1);
//   // CGAL::draw(arr);
// }

// void draw_spherical_arr() {
//   using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
//   using Direction_3 = Exact_kernel::Direction_3;
//   using Geodesic_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Exact_kernel>;
//   using Spherical_topo_traits = CGAL::Arr_spherical_topology_traits_2<Geodesic_traits>;
//   using Arrangement = CGAL::Arrangement_on_surface_2<Geodesic_traits, Spherical_topo_traits>;
//   using Point_2 = Geodesic_traits::Point_2;

//   Arrangement arr;
//   auto traits = arr.geometry_traits();
//   auto cst_pt = traits->construct_point_2_object();
//   auto cst_param = traits->parameter_space_in_x_2_object();

//   Point_2 p1 = cst_pt(Direction_3(1, 0, 0));
// }

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_2;
// void write_polyline(std::string filename, const std::vector<Point>& points) {
//   std::ofstream ofs_index("/Users/shep/codes/aos_2_js_helper/shapes.txt");
//   ofs_index << filename << std::endl;
//   std::ofstream ofs("/Users/shep/codes/aos_2_js_helper/" + filename);
//   for(const auto& pt : points) {
//     ofs << pt << "\n";
//   }
//   ofs << std::endl;
// }

int main() {
  // draw_segments_arr_2();
  draw_linear_arr_2();
  // test_zone();
  // draw_conic_arcs_arr();
  // draw_algebraic_arr();
  // draw_rational_arr();
  // draw_circle_segs_arr();
  return 0;
}