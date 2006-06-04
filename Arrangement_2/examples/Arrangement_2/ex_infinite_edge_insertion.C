//! \file examples/Arrangement_2/ex_infinite_edge_insertion.C
// Constructing an arrangement of unbounded linear objects using the simple
// edge-insertion functions.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>

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
typedef CGAL::Arr_naive_point_location<Arrangement_2> Naive_pl;

int main ()
{
  Arrangement_2      arr;
  Naive_pl           naive_pl (arr);

  Ray_2              ray1 (Point_2 (0, 1), Point_2 (0, 2));
  X_monotone_curve_2 cv1 (ray1);
  Line_2             line2 (Point_2 (0, 0), Point_2 (1, 1));
  X_monotone_curve_2 cv2 (line2);
  Ray_2              ray3 (Point_2 (0, 1), Point_2 (1, 2));
  X_monotone_curve_2 cv3 (ray3);
  Ray_2              ray4 (Point_2 (0, 1), Point_2 (-1, 2));
  X_monotone_curve_2 cv4 (ray4);

  Halfedge_handle    e1 = arr.insert_in_face_interior (cv1,
                                                       arr.unbounded_face());
  Vertex_handle      v1 = e1->source();
  Halfedge_handle    e2 = arr.insert_in_face_interior (cv2,
                                                       arr.unbounded_face());
  Halfedge_handle    e3 = arr.insert_from_left_vertex (cv3, v1);
  Halfedge_handle    e4 = arr.insert_from_right_vertex (cv4, v1);

  // Print out the size of the resulting arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << " (" << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  CGAL::Object                                 obj;
  Arrangement_2::Face_const_handle             f;
  Arrangement_2::Ccb_halfedge_const_circulator first, circ;

  obj = naive_pl.locate (Point_2 (1, 3));
  if (CGAL::assign (f ,obj))
  {
    std::cout << "Face is " << (f->is_unbounded() ? "unbounded." : "bounded.")
              << std::endl;
    first = circ = f->outer_ccb();
    do
    {
      if (! circ->is_fictitious())
        std::cout << "  " << circ->curve() << std::endl;
      ++circ;
    } while (circ != first);
  }

  obj = naive_pl.locate (Point_2 (3, -3));
  if (CGAL::assign (f ,obj))
  {
    std::cout << "Face is " << (f->is_unbounded() ? "unbounded." : "bounded.")
              << std::endl;
    first = circ = f->outer_ccb();
    do
    {
      if (! circ->is_fictitious())
        std::cout << "  " << circ->curve() << std::endl;
      ++circ;
    } while (circ != first);
  }

  return (0);
}
