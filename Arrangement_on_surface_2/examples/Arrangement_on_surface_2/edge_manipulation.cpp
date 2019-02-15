//! \file examples/Arrangement_on_surface_2/edge_manipulation.cpp
// Using the edge-manipulation functions.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include "arr_print.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                  Traits_2;
typedef Traits_2::Point_2                                   Point_2;
typedef Traits_2::X_monotone_curve_2                        Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                       Arrangement_2;
typedef Arrangement_2::Vertex_handle                        Vertex_handle;
typedef Arrangement_2::Halfedge_handle                      Halfedge_handle;

int main()
{
  // Step(a) - construct a triangular face.
  Arrangement_2   arr;

  Segment_2       s1(Point_2(667, 1000), Point_2(4000, 5000));
  Segment_2       s2(Point_2(4000, 0), Point_2(4000, 5000));
  Segment_2       s3(Point_2(667, 1000), Point_2(4000, 0));

  Halfedge_handle e1 = arr.insert_in_face_interior(s1, arr.unbounded_face());
  Vertex_handle   v1 = e1->source();
  Vertex_handle   v2 = e1->target();
  Halfedge_handle e2 = arr.insert_from_right_vertex(s2, v2);
  Vertex_handle   v3 = e2->target();
  arr.insert_at_vertices(s3, v3, v1);

  // Step (b) - create additional two faces inside the triangle.
  Point_2         p1(4000, 3666), p2(4000, 1000);
  Segment_2       s4(Point_2(4000, 5000), p1);
  Segment_2       s5(p1, p2);
  Segment_2       s6(Point_2(4000, 0), p2);

  Halfedge_handle e4 = arr.split_edge(e2, s4, Segment_2(Point_2(4000, 0), p1));
  Vertex_handle   w1 = e4->target();
  Halfedge_handle e5 = arr.split_edge(e4->next(), s5, s6);
  Vertex_handle   w2 = e5->target();
  Halfedge_handle e6 = e5->next();

  Segment_2       s7(p1, Point_2(3000, 2666));
  Segment_2       s8(p2, Point_2(3000, 1333));
  Segment_2       s9(Point_2(3000, 2666), Point_2(2000, 1666));
  Segment_2       s10(Point_2(3000, 1333), Point_2(2000, 1666));
  Segment_2       s11(Point_2(3000, 1333), Point_2(3000, 2666));

  Halfedge_handle e7 = arr.insert_from_right_vertex(s7, w1);
  Vertex_handle   v4 = e7->target();
  Halfedge_handle e8 = arr.insert_from_right_vertex(s8, w2);
  Vertex_handle   v5 = e8->target();
  Vertex_handle   v6 =
    arr.insert_in_face_interior(Point_2(2000, 1666), e8->face());

  arr.insert_at_vertices(s9, v4, v6);
  arr.insert_at_vertices(s10, v5, v6);
  arr.insert_at_vertices(s11, v4, v5);

  // Step(c) - remove and merge faces to form a single hole in the traingle.
  arr.remove_edge(e7);
  arr.remove_edge(e8);

  e5 = arr.merge_edge(e5, e6, Segment_2(e5->source()->point(),
                                        e6->target()->point()));
  e2 = arr.merge_edge(e4, e5, s2);

  print_arrangement(arr);
  return 0;
}
