//! \file examples/Arrangement_on_surface_2/ex_vertical_ray_shooting.cpp
// Answering vertical ray-shooting queries.

#include <CGAL/basic.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>

#include "arr_inexact_construction_segments.h"
#include "point_location_utils.h"

typedef CGAL::Arr_walk_along_line_point_location<Arrangement> Walk_pl;
typedef CGAL::Arr_trapezoid_ric_point_location<Arrangement>   Trap_pl;

int main() {
  // Construct the arrangement.
  Arrangement arr;
  construct_segments_arr(arr);

  // Perform some vertical ray-shooting queries using the walk strategy.
  Walk_pl walk_pl(arr);
  shoot_vertical_ray(walk_pl, Point(1, 4));
  shoot_vertical_ray(walk_pl, Point(4, 3));
  shoot_vertical_ray(walk_pl, Point(6, 3));

  // Attach the trapezoid-RIC object to the arrangement and perform queries.
  Trap_pl trap_pl(arr);
  shoot_vertical_ray(trap_pl, Point(3, 2));
  shoot_vertical_ray(trap_pl, Point(5, 2));
  shoot_vertical_ray(trap_pl, Point(1, 0));

  return 0;
}
