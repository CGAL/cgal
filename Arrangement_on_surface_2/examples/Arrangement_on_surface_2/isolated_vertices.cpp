//! \file examples/Arrangement_on_surface_2/isolated_vertices.cpp
// Constructing an arrangement with isolated vertices.

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
typedef Arrangement_2::Face_handle                    Face_handle;

int main ()
{
  Arrangement_2   arr;

  // Insert some isolated points:
  Face_handle     uf = arr.unbounded_face();
  Vertex_handle   u1 = arr.insert_in_face_interior (Point_2 (3, 3), uf);
  Vertex_handle   u2 = arr.insert_in_face_interior (Point_2 (1, 5), uf);
  Vertex_handle   u3 = arr.insert_in_face_interior (Point_2 (5, 5), uf);

  // Insert four segments that form a rectangular face:
  Segment_2       s1 (Point_2 (1, 3), Point_2 (3, 5));
  Segment_2       s2 (Point_2 (3, 5), Point_2 (5, 3));
  Segment_2       s3 (Point_2 (5, 3), Point_2 (3, 1));
  Segment_2       s4 (Point_2 (3, 1), Point_2 (1, 3));

  Halfedge_handle e1 = arr.insert_in_face_interior (s1, uf);
  Vertex_handle   v1 = e1->source();
  Vertex_handle   v2 = e1->target();
  Halfedge_handle e2 = arr.insert_from_left_vertex (s2, v2);
  Vertex_handle   v3 = e2->target();
  Halfedge_handle e3 = arr.insert_from_right_vertex (s3, v3);
  Vertex_handle   v4 = e3->target();
  Halfedge_handle e4 = arr.insert_at_vertices (s4, v4, v1);

  // Remove the isolated vertices located in the unbounded face.
  Arrangement_2::Vertex_iterator        curr_v, next_v;

  for (curr_v = arr.vertices_begin();
       curr_v != arr.vertices_end(); curr_v = next_v)
  {
    // Store an iterator to the next vertex (as we may delete curr_v and
    // invalidate the iterator).
    next_v = curr_v;
    ++next_v;

    if (curr_v->is_isolated() && curr_v->face() == uf)
      arr.remove_isolated_vertex (curr_v);      
  }

  print_arrangement (arr);
  return 0;
}
