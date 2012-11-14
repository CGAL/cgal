//! \file examples/Arrangement_on_surface_2/special_edge_insertion.cpp
// Constructing an arrangement using the specialized edge-insertion functions.

#include "arr_inexact_construction_segments.h"
#include "arr_print.h"

int main()
{
  Point           p0(3, 3), p1(1, 3), p2(3, 5), p3(5, 3), p4(3, 1);
  Segment         s1(p1, p2), s2(p2, p3), s3(p3, p4), s4(p4, p1);
  Segment         s5(p1, p0), s6(p0, p3), s7(p4, p0), s8(p0, p2);

  Arrangement     arr;
  Vertex_handle   v0 = arr.insert_in_face_interior(p0, arr.unbounded_face());
  Halfedge_handle e1 = arr.insert_in_face_interior(s1, arr.unbounded_face());
  Halfedge_handle e2 = arr.insert_from_left_vertex(s2, e1);
  Halfedge_handle e3 = arr.insert_from_right_vertex(s3, e2);
  Halfedge_handle e4 = arr.insert_at_vertices(s4, e3, e1->twin());
  Halfedge_handle e5 = arr.insert_at_vertices(s5, e1->twin(), v0);
  Halfedge_handle e6 = arr.insert_at_vertices(s6, e5, e3->twin());
  arr.insert_at_vertices(s7, e4->twin(), e6->twin());
  arr.insert_at_vertices(s8, e5, e2->twin());

  print_arrangement(arr);
  return 0;
}
