//! \file examples/Arrangement_on_surface_2/isolated_vertices.cpp
// Constructing an arrangement with isolated vertices.

#include "arr_inexact_construction_segments.h"
#include "arr_print.h"

int main()
{
  // Insert isolated points.
  Arrangement arr;
  Face_handle uf = arr.unbounded_face();
  arr.insert_in_face_interior(Point(3, 3), uf);
  arr.insert_in_face_interior(Point(1, 5), uf);
  arr.insert_in_face_interior(Point(5, 5), uf);

  // Insert four segments that form a square-shaped face.
  Point p1(1, 3), p2(3, 5), p3(5, 3), p4(3, 1);
  Segment s1(p1, p2), s2(p2, p3), s3(p3, p4), s4(p4, p1);

  Halfedge_handle e1 = arr.insert_in_face_interior(s1, uf);
  Vertex_handle   v1 = e1->source();
  Vertex_handle   v2 = e1->target();
  Halfedge_handle e2 = arr.insert_from_left_vertex(s2, v2);
  Vertex_handle   v3 = e2->target();
  Halfedge_handle e3 = arr.insert_from_right_vertex(s3, v3);
  Vertex_handle   v4 = e3->target();
  arr.insert_at_vertices(s4, v4, v1);

  // Remove the isolated vertices located in the unbounded face.
  Arrangement::Vertex_iterator curr, next = arr.vertices_begin();
  for (curr = next++; curr != arr.vertices_end(); curr = next++) {
    // Keep an iterator to the next vertex, as curr might be deleted.
    if (curr->is_isolated() && curr->face() == uf)
      arr.remove_isolated_vertex(curr);
  }

  print_arrangement(arr);
  return 0;
}
