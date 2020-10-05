//! \file examples/Arrangement_on_surface_2/point_location.cpp
// Answering point-location queries.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>

#include "point_location_utils.h"

typedef int                                                     Number_type;
typedef CGAL::Simple_cartesian<Number_type>                     Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Point_2                                       Point_2;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2>           Naive_pl;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2>       Landmarks_pl;

int main ()
{
  // Construct the arrangement.
  Arrangement_2    arr;
  Naive_pl         naive_pl(arr);
  construct_segments_arr(arr);

  // Perform some point-location queries using the naive strategy.
  point_location_query (naive_pl, Point_2(1, 4));        // q1
  point_location_query (naive_pl, Point_2(4, 3));        // q2
  point_location_query (naive_pl, Point_2(6, 3));        // q3

  // Attach the landmarks object to the arrangement and perform queries.
  Landmarks_pl landmarks_pl;
  landmarks_pl.attach(arr);
  point_location_query (landmarks_pl, Point_2(3, 2));    // q4
  point_location_query (landmarks_pl, Point_2(5, 2));    // q5
  point_location_query (landmarks_pl, Point_2(1, 0));    // q6

  return 0;
}
