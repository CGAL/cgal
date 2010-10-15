//! \file examples/Arrangement_on_surface_2/polylines.cpp
// Constructing an arrangement of polylines.

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>

#include "arr_print.h"

typedef CGAL::Quotient<CGAL::MP_Float>                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>   Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Polyline_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;

int main ()
{
  Arrangement_2         arr;

  Point_2               points1[5];
  points1[0] = Point_2 (0, 0);
  points1[1] = Point_2 (2, 4);
  points1[2] = Point_2 (3, 0);
  points1[3] = Point_2 (4, 4);
  points1[4] = Point_2 (6, 0);
  Polyline_2            pi1 (&points1[0], &points1[5]);

  std::list<Point_2>    points2;
  points2.push_back (Point_2 (1, 3));
  points2.push_back (Point_2 (0, 2));
  points2.push_back (Point_2 (1, 0));
  points2.push_back (Point_2 (2, 1));
  points2.push_back (Point_2 (3, 0));
  points2.push_back (Point_2 (4, 1));
  points2.push_back (Point_2 (5, 0));
  points2.push_back (Point_2 (6, 2));
  points2.push_back (Point_2 (5, 3));
  points2.push_back (Point_2 (4, 2));
  Polyline_2            pi2 (points2.begin(), points2.end());

  std::vector<Point_2>  points3 (4);
  points3[0] = Point_2 (0, 2);
  points3[1] = Point_2 (1, 2);
  points3[2] = Point_2 (3, 6);
  points3[3] = Point_2 (5, 2);
  Polyline_2            pi3 (points3.begin(), points3.end());
  
  insert (arr, pi1);
  insert (arr, pi2);
  insert (arr, pi3);
  
  print_arrangement (arr);
  return 0;
}
