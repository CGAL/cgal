#include <iostream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
// #include <CGAL/draw_arrangement_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Circle_2 = Kernel::Circle_2;

int main() {
  Kernel kernel;
  auto ctr_circle = kernel.construct_circle_2_object();
  auto circle1 = ctr_circle(Point_2(0, 1), 1);
  auto circle2 = ctr_circle(Point_2(0, -1), 1);
  auto circle3 = ctr_circle(Point_2(0, 2), 4);

  // 1. Circular arcs and linear segments
  using Circle_segment_arr_traits_2 = CGAL::Arr_circle_segment_traits_2<Kernel>;
  using Circle_segment_xcv_2 = Circle_segment_arr_traits_2::X_monotone_curve_2;
  using Circle_segment_pnt_2 = Circle_segment_arr_traits_2::Point_2;
  using Circle_segment_gps_traits_2 = CGAL::Gps_traits_2<Circle_segment_arr_traits_2>;
  using Circle_segment_polygon = Circle_segment_gps_traits_2::General_polygon_2;

  Circle_segment_arr_traits_2 circle_segment_traits;

  Circle_segment_pnt_2 cs_pnt11(1, 1);
  Circle_segment_pnt_2 cs_pnt12(-1, 1);
  Circle_segment_xcv_2 xcv11(circle1, cs_pnt11, cs_pnt12, CGAL::COUNTERCLOCKWISE);
  Circle_segment_xcv_2 xcv12(circle1, cs_pnt12, cs_pnt11, CGAL::COUNTERCLOCKWISE);
  Circle_segment_polygon pgn1;
  pgn1.push_back(xcv11);
  pgn1.push_back(xcv12);

  Circle_segment_pnt_2 cs_pnt21(1, -1);
  Circle_segment_pnt_2 cs_pnt22(-1, -1);
  Circle_segment_xcv_2 xcv21(circle2, cs_pnt21, cs_pnt22, CGAL::COUNTERCLOCKWISE);
  Circle_segment_xcv_2 xcv22(circle2, cs_pnt22, cs_pnt21, CGAL::COUNTERCLOCKWISE);
  Circle_segment_polygon pgn2;
  pgn2.push_back(xcv21);
  pgn2.push_back(xcv22);

  Circle_segment_pnt_2 cs_pnt31(2, 2);
  Circle_segment_pnt_2 cs_pnt32(-2, 2);
  Circle_segment_xcv_2 xcv31(circle3, cs_pnt31, cs_pnt32, CGAL::COUNTERCLOCKWISE);
  Circle_segment_xcv_2 xcv32(circle3, cs_pnt32, cs_pnt31, CGAL::COUNTERCLOCKWISE);
  Circle_segment_polygon pgn3;
  pgn3.push_back(xcv31);
  pgn3.push_back(xcv32);

  // 1.1.
  auto do_intersect = CGAL::do_intersect(pgn1, pgn2);
  if (do_intersect) {
    std::cerr << "The circles intersect (case 1)\n" << std::endl;
    return 1;
  }

  // 1.2.
  std::vector<Circle_segment_polygon> pgns1 = { pgn1, pgn2 };
  do_intersect = CGAL::do_intersect(pgns1.begin(), pgns1.end());
  if (do_intersect) {
    std::cerr << "The circles intersect (case 2)\n" << std::endl;
    return 1;
  }

  // 2.1.
  do_intersect = CGAL::do_intersect(pgn1, pgn3);
  if (! do_intersect) {
    std::cerr << "The circles do not intersect (case 1)\n" << std::endl;
    return 1;
  }

  // 2.2.
  std::vector<Circle_segment_polygon> pgns2 = { pgn1, pgn3 };
  do_intersect = CGAL::do_intersect(pgns2.begin(), pgns2.end());
  if (! do_intersect) {
    std::cerr << "The circles do not intersect (case 2)\n" << std::endl;
    return 1;
  }

  // using Circle_segment_arr = CGAL::Arrangement_2<Circle_segment_arr_traits_2>;
  // Circle_segment_arr arr;
  // CGAL::insert_non_intersecting_curve(arr, xcv11);
  // CGAL::insert_non_intersecting_curve(arr, xcv12);
  // CGAL::insert_non_intersecting_curve(arr, xcv21);
  // CGAL::insert_non_intersecting_curve(arr, xcv22);
  // CGAL::draw(arr);

  return 0;
}
