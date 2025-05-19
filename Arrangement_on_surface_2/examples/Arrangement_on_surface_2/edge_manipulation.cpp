//! \file examples/Arrangement_on_surface_2/edge_manipulation.cpp
// Using the edge-manipulation functions.

#include "arr_inexact_construction_segments.h"
#include "arr_print.h"

int main() {
  // Step (a)---construct a rectangular face.
  Point q1(1, 3), q2(3, 5), q3(5, 3), q4(3, 1);
  Segment s4(q1, q2), s1(q2, q3), s3(q3, q4), s2(q4, q1);

  Arrangement arr;
  Halfedge_handle e1 = arr.insert_in_face_interior(s1, arr.unbounded_face());
  Halfedge_handle e2 = arr.insert_in_face_interior(s2, arr.unbounded_face());

  e2 = e2->twin();     // as we wish e2 to be directed from right to left
  arr.insert_at_vertices(s3, e1->target(), e2->source());
  arr.insert_at_vertices(s4, e2->target(), e1->source());
  std::cout << "After step (a):\n";
  print_arrangement(arr);

  // Step (b)---split e1 and e2 and connect the split points with a segment.
  Point p1(4,4), p2(2,2);
  Segment s1_1(q2, p1), s1_2(p1, q3), s2_1(q4, p2), s2_2(p2, q1), s(p1, p2);

  e1 = arr.split_edge(e1, s1_1, s1_2);
  e2 = arr.split_edge(e2, s2_1, s2_2);
  Halfedge_handle e = arr.insert_at_vertices(s, e1->target(), e2->target());
  std::cout << std::endl << "After step (b):" << std::endl;
  print_arrangement(arr);

  // Step (c)---remove the edge e and merge e1 and e2 with their successors.
  arr.remove_edge(e);
  arr.merge_edge(e1, e1->next(), s1);
  arr.merge_edge(e2, e2->next(), s2);
  std::cout << std::endl << "After step (c):\n";
  print_arrangement(arr);
  return 0;
}
