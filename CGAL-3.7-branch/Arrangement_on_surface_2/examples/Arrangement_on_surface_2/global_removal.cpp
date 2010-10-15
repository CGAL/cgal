//! \file examples/Arrangement_on_surface_2/global_removal.cpp
// Using the global removal functions.

#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>

#include "arr_print.h"

typedef int                                           Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;
typedef CGAL::Arr_naive_point_location<Arrangement_2> Naive_pl;

int main ()
{
  // Create an arrangement of four line segments forming an H-shape:
  Arrangement_2   arr;
  Naive_pl        pl (arr);

  Segment_2       s1 (Point_2 (1, 3), Point_2 (4, 3));
  Halfedge_handle e1 = arr.insert_in_face_interior (s1, arr.unbounded_face()); 
  Segment_2       s2 (Point_2 (1, 4), Point_2 (4, 4));
  Halfedge_handle e2 = arr.insert_in_face_interior (s2, arr.unbounded_face()); 
  Segment_2       s3 (Point_2 (1, 1), Point_2 (1, 6));
  Segment_2       s4 (Point_2 (4, 1), Point_2 (4, 6));

  insert (arr, s3, pl);
  insert (arr, s4, pl);

  std::cout << "The initial arrangement:" << std::endl;
  print_arrangement (arr);

  // Remove the horizontal edge from the arrangement, and its end vertices:
  Vertex_handle   v1 = e1->source(), v2 = e1->target();
  arr.remove_edge (e1);
  remove_vertex (arr, v1);
  remove_vertex (arr, v2);

  // Remove the second horizontal edge e2 from the arrangement:
  remove_edge (arr, e2);

  std::cout << "The final arrangement:" << std::endl;
  print_arrangement (arr);
  return 0;
}
