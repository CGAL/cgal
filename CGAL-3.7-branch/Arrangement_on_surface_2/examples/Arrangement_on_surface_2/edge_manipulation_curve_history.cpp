//! \file examples/Arrangement_on_surface_2/edge_manipulation_curve_history.cpp
// Removing curves and manipulating edges in an arrangement with history.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_with_history_2.h>

typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef Kernel::Point_2                               Rat_point_2;
typedef Kernel::Circle_2                              Circle_2;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>     Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Curve_2;
typedef CGAL::Arrangement_with_history_2<Traits_2>    Arr_with_hist_2;
typedef Arr_with_hist_2::Curve_handle                 Curve_handle;
typedef CGAL::Arr_walk_along_line_point_location<Arr_with_hist_2>  
                                                      Point_location;

int main ()
{
  // Construct an arrangement containing nine circles: C[0] of radius 2 and
  // C[1], ..., C[8] of radius 1.
  const Number_type _7_halves = Number_type (7, 2); 
  Arr_with_hist_2   arr;
  Curve_2           C[9];
  Curve_handle      handles[9];
  int               k;

  C[0] = Circle_2 (Rat_point_2 (_7_halves, _7_halves), 4, CGAL::CLOCKWISE);
  C[1] = Circle_2 (Rat_point_2 (_7_halves, 6), 1, CGAL::CLOCKWISE);
  C[2] = Circle_2 (Rat_point_2 (5, 6), 1, CGAL::CLOCKWISE);
  C[3] = Circle_2 (Rat_point_2 (6, _7_halves), 1, CGAL::CLOCKWISE);
  C[4] = Circle_2 (Rat_point_2 (5, 2), 1, CGAL::CLOCKWISE);
  C[5] = Circle_2 (Rat_point_2 (_7_halves, 1), 1, CGAL::CLOCKWISE);
  C[6] = Circle_2 (Rat_point_2 (2, 2), 1, CGAL::CLOCKWISE);
  C[7] = Circle_2 (Rat_point_2 (1, _7_halves), 1, CGAL::CLOCKWISE);
  C[8] = Circle_2 (Rat_point_2 (2, 5), 1, CGAL::CLOCKWISE);

  for (k = 0; k < 9; k++)
    handles[k] = insert (arr, C[k]);

  std::cout << "The initial arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Remove the large circle C[0].
  std::cout << "Removing C[0] : ";
  std::cout << remove_curve (arr, handles[0]) 
            << " edges have been removed." << std::endl;

  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Locate the point q, which should be on an edge e.
  Point_location                          pl (arr);
  const Point_2                           q = Point_2 (_7_halves, 7);
  CGAL::Object                            obj = pl.locate (q);
  Arr_with_hist_2::Halfedge_const_handle  e;

  CGAL_assertion_code(bool success = ) CGAL::assign (e, obj);
  CGAL_assertion (success);
 
  // Split the edge e to two edges e1 and e2;
  Arr_with_hist_2::Halfedge_handle        e1, e2;

  e1 = arr.split_edge (arr.non_const_handle (e), q);
  e2 = e1->next();

  std::cout << "After edge split: "
            << "V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Merge back the two split edges.
  arr.merge_edge (e1, e2);

  std::cout << "After edge merge: "
            << "V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;
  return 0;
}
