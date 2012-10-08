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
  Naive_pl         naive_pl (arr);
  Landmarks_pl     landmarks_pl;

  construct_segments_arr (arr);

  // Perform some point-location queries using the naive strategy.
  Point_2          q1 (1, 4);
  Point_2          q2 (4, 3);
  Point_2          q3 (6, 3);

  point_location_query (naive_pl, q1);
  point_location_query (naive_pl, q2);
  point_location_query (naive_pl, q3);

  // Attach the landmarks object to the arrangement and perform queries.
  Point_2          q4 (3, 2);
  Point_2          q5 (5, 2);
  Point_2          q6 (1, 0);

  landmarks_pl.attach (arr);

  point_location_query (landmarks_pl, q4);
  point_location_query (landmarks_pl, q5);
  point_location_query (landmarks_pl, q6);
  
  return 0;
}
