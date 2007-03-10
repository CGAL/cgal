//! \file examples/Arrangement_2/ex_unbounded_non_intersecting.cpp
// Constructing an arrangement of unbounded linear objects using the insertion
// function for non-intersecting curves.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef int                                           Number_type;
typedef CGAL::Simple_cartesian<Number_type>           Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>             Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Segment_2                           Segment_2;
typedef Traits_2::Ray_2                               Ray_2;
typedef Traits_2::Line_2                              Line_2;
typedef Traits_2::X_monotone_curve_2                  X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;

int main ()
{
  Arrangement_2      arr;

  // Insert a line in the (currently single) unbounded face of the arrangement,
  // then split it into two at (0,0). Assign v to be the split point.
  X_monotone_curve_2 c1 = Line_2 (Point_2 (-1, 0), Point_2 (1, 0));
  Halfedge_handle    e1 = arr.insert_in_face_interior (c1,
                                                       arr.unbounded_face());

  X_monotone_curve_2 c1_left = Ray_2 (Point_2 (0, 0), Point_2 (-1, 0));
  X_monotone_curve_2 c1_right = Ray_2 (Point_2 (0, 0), Point_2 (1, 0));

  e1 = arr.split_edge (e1, c1_left, c1_right);
  Vertex_handle      v = e1->target();

  CGAL_assertion (! v->is_at_infinity());

  // Add two more rays using the specialized insertion functions.
  X_monotone_curve_2 c2 = Ray_2 (Point_2 (0, 0), Point_2 (-1, 1));
  X_monotone_curve_2 c3 = Ray_2 (Point_2 (0, 0), Point_2 (1, 1));

  arr.insert_from_right_vertex (c2, v);
  arr.insert_from_left_vertex (c3, v);

  // Insert three more interior-disjoint rays.
  X_monotone_curve_2 c4 = Ray_2 (Point_2 (0, -1), Point_2 (-2, -2));
  X_monotone_curve_2 c5 = Ray_2 (Point_2 (0, -1), Point_2 (2, -2));
  X_monotone_curve_2 c6 = Ray_2 (Point_2 (0, 0), Point_2 (0, 1));

  insert_non_intersecting_curve (arr, c4);
  insert_non_intersecting_curve (arr, c5);
  insert_non_intersecting_curve (arr, c6);

  // Print out the size of the resulting arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << " (plus " << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces()
            << " (" << arr.number_of_unbounded_faces() << " unbounded)"
            << std::endl << std::endl;

  return (0);
}
