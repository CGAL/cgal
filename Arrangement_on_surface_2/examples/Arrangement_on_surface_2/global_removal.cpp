//! \file examples/Arrangement_on_surface_2/global_removal.cpp
// Using the global removal functions.

#include "arr_exact_construction_segments.h"
#include "arr_print.h"

int main() {
  // Create an arrangement of four line segments forming an H-shape:
  Arrangement arr;

  Segment s1(Point(1, 3), Point(5, 3));
  Halfedge_handle e1 = arr.insert_in_face_interior(s1, arr.unbounded_face());
  Segment s2(Point(1, 4), Point(5, 4));
  Halfedge_handle e2 = arr.insert_in_face_interior (s2, arr.unbounded_face());
  insert(arr, Segment(Point(1, 1), Point(1, 6)));
  insert(arr, Segment(Point(5, 1), Point(5, 6)));

  std::cout << "The initial arrangement:\n";
  print_arrangement(arr);

  // Remove e1 and its incident vertices using the function remove_edge().
  Vertex_handle v1 = e1->source(), v2 = e1->target();
  arr.remove_edge(e1);
  remove_vertex(arr, v1);
  remove_vertex(arr, v2);

  // Remove e2 using the free remove_edge() function.
  remove_edge(arr, e2);

  std::cout << "The final arrangement:\n";
  print_arrangement(arr);
  return 0;
}
