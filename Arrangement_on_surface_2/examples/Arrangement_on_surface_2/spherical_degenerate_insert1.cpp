//! \file examples/Arrangement_on_surface_2/aggregated_insertion.cpp
// Using the global aggregated insertion functions.

#include "arr_geodesic_on_sphere.h"
#include "arr_print.h"

int main() {
  Geom_traits traits;
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_xcv = traits.construct_x_monotone_curve_2_object();

  Arrangement arr(&traits);
  Point sp = ctr_p(0, 0, -1);
  Point np = ctr_p(0, 0, 1);
  Vertex_handle spv = arr.insert_in_face_interior(sp, arr.reference_face());
  Vertex_handle npv = arr.insert_in_face_interior(np, arr.reference_face());

#if 1
  Point p1 = ctr_p(-1,  0, 0);
#else
  Point p1 = ctr_p(-1,  -1, 0);
#endif
  Point p2 = ctr_p( 0, -1, 0);
  Point p3 = ctr_p( 1,  0, 0);

  Vertex_handle v1 = arr.insert_in_face_interior(p1, arr.reference_face());
  Vertex_handle v2 = arr.insert_in_face_interior(p2, arr.reference_face());
  Vertex_handle v3 = arr.insert_in_face_interior(p3, arr.reference_face());

  X_monotone_curve xcv_sp1 = ctr_xcv(sp, p1);
  X_monotone_curve xcv_sp2 = ctr_xcv(sp, p2);
  X_monotone_curve xcv_sp3 = ctr_xcv(sp, p3);
  arr.insert_at_vertices(xcv_sp1, spv, v1);
  arr.insert_at_vertices(xcv_sp2, spv, v2);
  arr.insert_at_vertices(xcv_sp3, spv, v3);

  X_monotone_curve xcv_np1 = ctr_xcv(np, p1);
  X_monotone_curve xcv_np2 = ctr_xcv(np, p2);
  X_monotone_curve xcv_np3 = ctr_xcv(np, p3);
  arr.insert_at_vertices(xcv_np1, npv, v1);
  arr.insert_at_vertices(xcv_np2, npv, v2);
  arr.insert_at_vertices(xcv_np3, npv, v3);

  X_monotone_curve xcv_12 = ctr_xcv(p1, p2);
  X_monotone_curve xcv_23 = ctr_xcv(p2, p3);
  arr.insert_at_vertices(xcv_12, v1, v2);
  arr.insert_at_vertices(xcv_23, v2, v3);

  print_arrangement_size(arr);          // print the arrangement size

  return 0;
}
