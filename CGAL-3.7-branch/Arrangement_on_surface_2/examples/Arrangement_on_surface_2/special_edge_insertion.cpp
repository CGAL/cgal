//! \file examples/Arrangement_on_surface_2/special_edge_insertion.cpp
// Constructing an arrangement using the specialized edge-insertion functions.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include "arr_print.h"

typedef int                                           Number_type;
typedef CGAL::Simple_cartesian<Number_type>           Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;

int main ()
{
  Arrangement_2   arr;

  Point_2         p0 (3, 3);
  Point_2         p1 (1, 3), p2 (3, 5), p3 (5, 3), p4 (3, 1);
  Segment_2       s1 (p1, p2);
  Segment_2       s2 (p2, p3);
  Segment_2       s3 (p3, p4);
  Segment_2       s4 (p4, p1);
  Segment_2       s5 (p1, p0);
  Segment_2       s6 (p0, p3);
  Segment_2       s7 (p4, p0);
  Segment_2       s8 (p0, p2);

  Vertex_handle   v0 = arr.insert_in_face_interior (p0, arr.unbounded_face());
  Halfedge_handle e1 = arr.insert_in_face_interior (s1, arr.unbounded_face());
  Halfedge_handle e2 = arr.insert_from_left_vertex (s2, e1);
  Halfedge_handle e3 = arr.insert_from_right_vertex (s3, e2);
  Halfedge_handle e4 = arr.insert_at_vertices (s4, e3, e1->twin());
  Halfedge_handle e5 = arr.insert_at_vertices (s5, e1->twin(), v0);
  Halfedge_handle e6 = arr.insert_at_vertices (s6, e5, e3->twin());
  Halfedge_handle e7 = arr.insert_at_vertices (s7, e4->twin(), e6->twin());
  Halfedge_handle e8 = arr.insert_at_vertices (s8, e5, e2->twin());

  print_arrangement (arr);
  return 0;
}
