//! \file examples/Arrangement_on_surface_2/point_location.cpp
// Answering point-location queries.

#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>

#include "arr_inexact_construction_segments.h"
#include "point_location_utils.h"

using Naive_pl = CGAL::Arr_naive_point_location<Arrangement>;
using Landmarks_pl = CGAL::Arr_landmarks_point_location<Arrangement>;

int main() {
  // Construct the arrangement.
  Arrangement arr;
  construct_segments_arr(arr);

  // Perform some point-location queries using the naive strategy.
  Naive_pl naive_pl(arr);
  locate_point(naive_pl, Point(1, 4));          // q1
  locate_point(naive_pl, Point(4, 3));          // q2
  locate_point(naive_pl, Point(6, 3));          // q3

  // Perform some point-location queries using the landmark strategy.
  Landmarks_pl landmarks_pl(arr);
  locate_point(landmarks_pl, Point(3, 2));      // q4
  locate_point(landmarks_pl, Point(5, 2));      // q5
  locate_point(landmarks_pl, Point(1, 0));      // q6

  return 0;
}
