//! \file examples/Arrangement_on_surface_2/incremental_insertion.cpp
// Using the global incremental insertion functions.

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

#include "arr_print.h"

typedef CGAL::Quotient<int>                                     Number_type;
typedef CGAL::Cartesian<Number_type>                            Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Point_2                                       Point_2;
typedef Traits_2::X_monotone_curve_2                            Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;

int main()
{
  // Construct the arrangement of five intersecting segments.
  Arrangement_2  arr;
  Walk_pl        pl(arr);

  Segment_2      s1(Point_2(1, 0), Point_2(2, 4));
  Segment_2      s2(Point_2(5, 0), Point_2(5, 5));
  Segment_2      s3(Point_2(1, 0), Point_2(5, 3));
  Segment_2      s4(Point_2(0, 2), Point_2(6, 0));
  Segment_2      s5(Point_2(3, 0), Point_2(5, 5));

  insert_non_intersecting_curve(arr, s1, pl);
  insert_non_intersecting_curve(arr, s2, pl);
  insert(arr, s3, pl);
  insert(arr, s4, pl);
  insert(arr, s5, pl);

  // Print the size of the arrangement.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Perform a point-location query on the resulting arrangement and print
  // the boundary of the face that contains it.
  Point_2 q(4, 1);
  Walk_pl::result_type obj = pl.locate(q);

  Arrangement_2::Face_const_handle  f;
  CGAL_assertion_code(bool success =) CGAL::assign(f, obj);

  CGAL_assertion(success);
  std::cout << "The query point (" << q << ") is located in: ";
  print_face<Arrangement_2>(f);

  return 0;
}
