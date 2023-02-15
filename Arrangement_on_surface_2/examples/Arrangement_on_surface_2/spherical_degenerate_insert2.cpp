//! \file examples/Arrangement_on_surface_2/aggregated_insertion.cpp
// Using the global aggregated insertion functions.

// #define CGAL_IDENTIFICATION_XY CGAL_X_MINUS_11_Y_7
// #define CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE 1

#include "arr_geodesic_on_sphere.h"
#include "arr_print.h"

int main() {
  Geom_traits traits;
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_xcv = traits.construct_x_monotone_curve_2_object();
  Arrangement arr(&traits);
  Point sp = ctr_p(0, 0, -1);
  Point np = ctr_p(0, 0, 1);

  Point p1 = ctr_p(-1,  0, -1);
  Point p2 = ctr_p(-1,  0,  1);
  X_monotone_curve xcv_sp_p2 = ctr_xcv(sp, p2);
  X_monotone_curve xcv_np_p1 = ctr_xcv(np, p1);
  // std::cout << "Inserting " << xcv_sp_p2 << std::endl;
  insert(arr, xcv_sp_p2);
  // std::cout << "Inserting " << xcv_np_p1 << std::endl;
  insert(arr, xcv_np_p1);

  Point q1 = ctr_p(-1,  -1, -1);
  Point q2 = ctr_p(-1,  -1,  1);
  X_monotone_curve xcv_sp_q2 = ctr_xcv(sp, q2);
  X_monotone_curve xcv_np_q1 = ctr_xcv(np, q1);
  insert(arr, xcv_sp_q2);
  insert(arr, xcv_np_q1);

  X_monotone_curve xcv_p1_q1 = ctr_xcv(p1, q1);
  X_monotone_curve xcv_p2_q2 = ctr_xcv(p2, q2);
  insert(arr, xcv_p1_q1);
  insert(arr, xcv_p2_q2);

  print_arrangement_size(arr);          // print the arrangement size

  return 0;
}
