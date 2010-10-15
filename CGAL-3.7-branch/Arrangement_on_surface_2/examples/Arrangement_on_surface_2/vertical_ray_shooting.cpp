//! \file examples/Arrangement_on_surface_2/ex_vertical_ray_shooting.cpp
// Answering vertical ray-shooting queries.

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>

#include "point_location_utils.h"

typedef CGAL::MP_Float                                          Number_type;
typedef CGAL::Cartesian<Number_type>                            Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Point_2                                       Point_2;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;
typedef CGAL::Arr_trapezoid_ric_point_location<Arrangement_2>   Trap_pl;

int main ()
{
  // Construct the arrangement.
  Arrangement_2    arr;
  Walk_pl          walk_pl (arr);
  Trap_pl          trap_pl;

  construct_segments_arr (arr);

  // Perform some vertical ray-shooting queries using the walk strategy.
  Point_2          q1 (1, 4);
  Point_2          q2 (4, 3);
  Point_2          q3 (6, 3);

  vertical_ray_shooting_query (walk_pl, q1);
  vertical_ray_shooting_query (walk_pl, q2);
  vertical_ray_shooting_query (walk_pl, q3);

  // Attach the trapezoid-RIC object to the arrangement and perform queries.
  Point_2          q4 (3, 2);
  Point_2          q5 (5, 2);
  Point_2          q6 (1, 0);

  trap_pl.attach (arr);
  vertical_ray_shooting_query (trap_pl, q4);
  vertical_ray_shooting_query (trap_pl, q5);
  vertical_ray_shooting_query (trap_pl, q6);

  return 0;
}
