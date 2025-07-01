
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
#include <vector>

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

// void draw_segments_arr_4() {
//   using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
//   using Segment_traits = CGAL::Arr_segment_traits_2<Exact_kernel>;
//   using Point_2 = Segment_traits::Point_2;
//   using X_monotone_curve_2 = Segment_traits::X_monotone_curve_2;
//   using Arrangement = CGAL::Arrangement_2<Segment_traits>;
//   Arrangement arr;
//   auto traits = arr.traits();
//   auto cst_x_curve = traits->construct_x_monotone_curve_2_object();

//   std::vector<X_monotone_curve_2> segments;

//   std::vector<Point_2> polyline1{
//       {6, -21}, {6, -3}, {13, -3}, {13, -12}, {18, -12}, {18, -4}, {18, -3}, {25, -3}, {25, -21}, {6, -21},
//   };
//   for(size_t i = 0; i < polyline1.size() - 1; ++i) {
//     segments.push_back(cst_x_curve(polyline1[i], polyline1[i + 1]));
//   }

//   std::vector<Point_2> polyline2{
//       {-27, -14}, {-24, -14}, {-21, -16}, {-19, -18}, {-18, -21}, {-17, -24}, {-18, -28}, {-19, -30}, {-22, -32},
//       {-24, -33}, {-29, -33}, {-34, -33}, {-38, -32}, {-43, -30}, {-46, -29}, {-50, -25}, {-53, -21}, {-53, -17},
//       {-54, -12}, {-53, -7},  {-52, -2},  {-50, 2},   {-45, 5},   {-40, 6},   {-34, 7},   {-29, 7},   {-23, 7},
//       {-20, 6},   {-18, 5},   {-16, 2},   {-15, 0},   {-16, -3},  {-17, -3},  {-18, -1},  {-19, 1},   {-20, 4},
//       {-22, 4},   {-26, 4},   {-28, 4},   {-31, 4},   {-33, 4},   {-37, 4},   {-41, 3},   {-44, 2},   {-47, 1},
//       {-48, 0},   {-49, -2},  {-50, -5},  {-51, -8},  {-51, -13}, {-51, -16}, {-50, -18}, {-49, -22}, {-46, -24},
//       {-42, -27}, {-40, -29}, {-36, -29}, {-32, -30}, {-28, -31}, {-24, -31}, {-21, -30}, {-19, -28}, {-19, -25},
//       {-20, -23}, {-21, -21}, {-24, -18}, {-26, -16}, {-27, -15}, {-27, -14},
//   };
//   for(size_t i = 0; i < polyline2.size() - 1; ++i) {
//     segments.push_back(cst_x_curve(polyline2[i], polyline2[i + 1]));
//   }

//   CGAL::insert(arr, segments.begin(), segments.end());

//   CGAL::draw_viewer(arr);
// }

// void draw_segments_arr_5() {
//   using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
//   using Segment_traits = CGAL::Arr_segment_traits_2<Exact_kernel>;
//   using Point_2 = Segment_traits::Point_2;
//   using X_monotone_curve_2 = Segment_traits::X_monotone_curve_2;
//   using Arrangement = CGAL::Arrangement_2<Segment_traits>;
//   Arrangement arr;
//   auto traits = arr.traits();
//   auto cst_x_curve = traits->construct_x_monotone_curve_2_object();

//   std::vector<X_monotone_curve_2> segments;

//   std::vector<Point_2> polyline1{
//       {11, 23},  {12, -13}, {17, -13}, {18, 1},  {18, 16},  {18, 23},  {15, 35},  {15, 44},  {18, 47},  {19, 42},
//       {19, 37},  {20, 28},  {23, 21},  {23, 9},  {24, -7},  {22, -15}, {20, -19}, {15, -19}, {10, -16}, {7, -14},
//       {7, -3},   {5, 11},   {1, 20},   {0, 35},  {1, 45},   {3, 51},   {5, 53},   {11, 53},  {19, 53},  {23, 49},
//       {25, 40},  {28, 31},  {34, 16},  {35, 7},  {38, -11}, {44, -11}, {43, -7},  {43, -5},  {42, 7},   {40, 23},
//       {38, 30},  {35, 49},  {19, 57},  {8, 58},  {-2, 55},  {-2, 30},  {-3, 15},  {-3, -16}, {1, -20},  {12, -22},
//       {20, -22}, {25, -21}, {26, -14}, {28, -3}, {27, 15},  {23, 31},  {22, 44},  {20, 49},  {15, 49},  {10, 46},
//       {8, 39},   {8, 32},   {8, 28},   {11, 23},
//   };
//   for(size_t i = 0; i < polyline1.size() - 1; ++i) {
//     segments.push_back(cst_x_curve(polyline1[i], polyline1[i + 1]));
//   }

//   CGAL::insert(arr, segments.begin(), segments.end());

//   CGAL::draw_viewer(arr);
// }

// void draw_segments_arr_6() {
//   using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
//   using Segment_traits = CGAL::Arr_segment_traits_2<Exact_kernel>;
//   using Point_2 = Segment_traits::Point_2;
//   using X_monotone_curve_2 = Segment_traits::X_monotone_curve_2;
//   using Arrangement = CGAL::Arrangement_2<Segment_traits>;
//   Arrangement arr;
//   auto traits = arr.traits();
//   auto cst_x_curve = traits->construct_x_monotone_curve_2_object();

//   std::vector<X_monotone_curve_2> segments;

//   std::vector<Point_2> polyline1{
//       {9, -3},   {11, -34}, {26, -34}, {26, -8}, {26, -34}, {36, -34}, {36, -31}, {36, -27}, {37, -23}, {41, -23},
//       {47, -23}, {49, -16}, {49, -2},  {42, 10}, {41, 0},   {41, -19}, {41, 0},   {42, 10},  {30, 10},  {28, 6},
//       {28, 2},   {29, -30}, {28, 2},   {28, 6},  {30, 10},  {42, 10},  {29, 16},  {18, 12},  {9, -3},
//   };
//   for(size_t i = 0; i < polyline1.size() - 1; ++i) {
//     segments.push_back(cst_x_curve(polyline1[i], polyline1[i + 1]));
//   }

//   CGAL::insert(arr, segments.begin(), segments.end());

//   CGAL::draw_viewer(arr);
// }

// void draw_linear_arr_1() {
//   using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
//   using Traits = CGAL::Arr_linear_traits_2<Exact_kernel>;
//   using Point_2 = Traits::Point_2;
//   using Line_2 = Traits::Line_2;
//   using Ray_2 = Traits::Ray_2;
//   using Curve_2 = Traits::Curve_2;
//   using Arrangement = CGAL::Arrangement_2<Traits>;
//   using Face_const_handle = Arrangement::Face_const_handle;
//   using Halfedge_const_handle = Arrangement::Halfedge_const_iterator;
//   using X_monotone_curve_2 = Traits::X_monotone_curve_2;

//   Arrangement arr;
//   auto x_axis = X_monotone_curve_2(Ray_2(Point_2(0, 0), Point_2(1, 0)));
//   auto y_axis = X_monotone_curve_2(Ray_2(Point_2(0, 0), Point_2(0, 1)));
//   CGAL::insert(arr, x_axis);
//   CGAL::insert(arr, y_axis);

//   CGAL::draw_viewer(arr);
// }

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
  int n = 5;
  for(int i = 0; i < n; ++i) {
    Point_2 p1(i * 5, 0);
    Point_2 p2(i * 5, 1);
    CGAL::insert(arr, Curve_2(Line_2(p1, p2)));
  }
  for(int i = 0; i < n; ++i) {
    Point_2 p1(0, i * 5);
    Point_2 p2(1, i * 5);
    CGAL::insert(arr, Curve_2(Line_2(p1, p2)));
  }
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
  std::cout << "Arrangement has " << arr.number_of_faces() << " faces." << std::endl;
  CGAL::draw_viewer(arr);
}

// void draw_linear_arr_3() {
//   using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
//   using Traits = CGAL::Arr_linear_traits_2<Exact_kernel>;
//   using Point_2 = Traits::Point_2;
//   using Line_2 = Traits::Line_2;
//   using Segment_2 = Traits::Segment_2;
//   using Ray_2 = Traits::Ray_2;
//   using Curve_2 = Traits::Curve_2;
//   using Arrangement = CGAL::Arrangement_2<Traits>;
//   using Face_const_handle = Arrangement::Face_const_handle;
//   using Halfedge_const_handle = Arrangement::Halfedge_const_iterator;
//   using X_monotone_curve_2 = Traits::X_monotone_curve_2;

//   std::vector<Point_2> points{
//       {11, 23},  {12, -13}, {17, -13}, {18, 1},  {18, 16},  {18, 23},  {15, 35},  {15, 44},  {18, 47},  {19, 42},
//       {19, 37},  {20, 28},  {23, 21},  {23, 9},  {24, -7},  {22, -15}, {20, -19}, {15, -19}, {10, -16}, {7, -14},
//       {7, -3},   {5, 11},   {1, 20},   {0, 35},  {1, 45},   {3, 51},   {5, 53},   {11, 53},  {19, 53},  {23, 49},
//       {25, 40},  {28, 31},  {34, 16},  {35, 7},  {38, -11}, {44, -11}, {43, -7},  {43, -5},  {42, 7},   {40, 23},
//       {38, 30},  {35, 49},  {19, 57},  {8, 58},  {-2, 55},  {-2, 30},  {-3, 15},  {-3, -16}, {1, -20},  {12, -22},
//       {20, -22}, {25, -21}, {26, -14}, {28, -3}, {27, 15},  {23, 31},  {22, 44},  {20, 49},  {15, 49},  {10, 46},
//       {8, 39},   {8, 32},   {8, 28},   {11, 23},
//   };

//   Arrangement arr;
//   auto& traits = *arr.traits();
//   std::vector<X_monotone_curve_2> segments;
//   for(size_t i = 0; i < points.size() - 1; ++i) {
//     Point_2 p1 = points[i];
//     Point_2 p2 = points[i + 1];
//     // create a segment
//     X_monotone_curve_2 seg = traits.construct_x_monotone_curve_2_object()(p1, p2);
//     segments.push_back(seg);
//   }

//   // insert segments into the arrangement
//   CGAL::insert(arr, segments.begin(), segments.end());

//   CGAL::draw_viewer(arr);
// }
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
  // draw_segments_arr_6();
  draw_linear_arr_2();
  // test_zone();
  // draw_conic_arcs_arr();
  // draw_algebraic_arr();
  // draw_rational_arr();
  // draw_circle_segs_arr();
  return 0;
}