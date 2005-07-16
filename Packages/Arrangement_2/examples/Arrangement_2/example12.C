// file: examples/Arrangement_2/example12.C

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include "arr_print.h"

typedef CGAL::MP_Float                                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>   Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::X_monotone_curve_2                    Polyline_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;

int main ()
{
  Arrangement_2 arr;

  Point_2 points1[5];
  points1[0] = Point_2(0, 0);
  points1[1] = Point_2(4, 8);
  points1[2] = Point_2(6, 0);
  points1[3] = Point_2(8, 8);
  points1[4] = Point_2(12, 0);
  Polyline_2 polyline1(&points1[0], &points1[5]);

  Point_2 points2[5];
  points2[0] = Point_2(0, 4);
  points2[1] = Point_2(2, 4);
  points2[2] = Point_2(6, 12);
  points2[3] = Point_2(10, 4);
  points2[4] = Point_2(12, 4);
  Polyline_2 polyline2(&points2[0], &points2[5]);

  Point_2 points3[11];
  points3[0] = Point_2(4, 4);
  points3[1] = Point_2(2, 6);
  points3[2] = Point_2(0, 4);
  points3[3] = Point_2(2, 0);
  points3[4] = Point_2(4, 2);
  points3[5] = Point_2(6, 0);
  points3[6] = Point_2(8, 2);
  points3[7] = Point_2(10, 0);
  points3[8] = Point_2(12, 4);
  points3[9] = Point_2(10, 6);
  points3[10] = Point_2(8, 4);
  Polyline_2 polyline3(&points3[0], &points3[11]);
  
  insert(arr, polyline1);
  insert(arr, polyline2);
  insert(arr, polyline3);
  
  print_arrangement (arr);
  return 0;
}
