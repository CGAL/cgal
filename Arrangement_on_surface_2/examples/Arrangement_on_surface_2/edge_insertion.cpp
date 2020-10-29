//! \file examples/Arrangement_on_surface_2/edge_insertion.cpp
// Constructing an arrangement using the simple edge-insertion functions.

#include "arr_inexact_construction_segments.h"
#include "arr_print.h"

int main() {
  Point p1(1, 3), p2(3, 5), p3(5, 3), p4(3, 1);
  Segment s1(p1, p2), s2(p2, p3), s3(p3, p4), s4(p4, p1), s5(p1, p3);

  Arrangement arr;
  Halfedge_handle e1 = arr.insert_in_face_interior(s1, arr.unbounded_face());
  Vertex_handle v1 = e1->source();
  Vertex_handle v2 = e1->target();
  Halfedge_handle e2 = arr.insert_from_left_vertex(s2, v2);
  Vertex_handle v3 = e2->target();
  Halfedge_handle e3 = arr.insert_from_right_vertex(s3, v3);
  Vertex_handle v4 = e3->target();
  arr.insert_at_vertices(s4, v4, v1);   // return e4
  arr.insert_at_vertices(s5, v1, v3);   // return e5

  print_arrangement(arr);
  return 0;
}
