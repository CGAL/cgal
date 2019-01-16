//! \file examples/Arrangement_on_surface_2/edge_insertion.cpp
// Constructing an arrangement using the simple edge-insertion functions.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include "arr_print.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT                                          Number_type;

typedef CGAL::Arr_segment_traits_2<Kernel>                  Traits;
typedef Traits::Point_2                                     Point;
typedef Traits::X_monotone_curve_2                          Segment;
typedef CGAL::Arrangement_2<Traits>                         Arrangement;
typedef Arrangement::Vertex_handle                          Vertex_handle;
typedef Arrangement::Halfedge_handle                        Halfedge_handle;

int main()
{
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
  arr.insert_at_vertices(s4, v4, v1);
  arr.insert_at_vertices(s5, v1, v3);

  print_arrangement(arr);
  return 0;
}
