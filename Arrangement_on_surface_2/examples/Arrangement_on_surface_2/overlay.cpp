//! \file examples/Arrangement_on_surface_2/overlay.cpp
// A simple overlay of two arrangements.

#include <CGAL/basic.h>
#include <CGAL/Arr_overlay_2.h>

#include "arr_exact_construction_segments.h"
#include "arr_print.h"

int main() {
  // Construct the first arrangement, containing a square-shaped face.
  Arrangement arr1;
  insert_non_intersecting_curve(arr1, Segment(Point(2, 2), Point(6, 2)));
  insert_non_intersecting_curve(arr1, Segment(Point(6, 2), Point(6, 6)));
  insert_non_intersecting_curve(arr1, Segment(Point(6, 6), Point(2, 6)));
  insert_non_intersecting_curve(arr1, Segment(Point(2, 6), Point(2, 2)));

  // Construct the second arrangement, containing a rhombus-shaped face.
  Arrangement arr2;
  insert_non_intersecting_curve(arr2, Segment(Point(4, 1), Point(7, 4)));
  insert_non_intersecting_curve(arr2, Segment(Point(7, 4), Point(4, 7)));
  insert_non_intersecting_curve(arr2, Segment(Point(4, 7), Point(1, 4)));
  insert_non_intersecting_curve(arr2, Segment(Point(1, 4), Point(4, 1)));

  // Compute the overlay of the two arrangements.
  Arrangement overlay_arr;
  CGAL::overlay(arr1, arr2, overlay_arr);
  print_arrangement_size(overlay_arr);

  return 0;
}
