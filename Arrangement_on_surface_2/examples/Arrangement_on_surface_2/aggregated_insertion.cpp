//! \file examples/Arrangement_on_surface_2/aggregated_insertion.cpp
// Using the global aggregated insertion functions.

#include "arr_exact_construction_segments.h"
#include "arr_print.h"

int main() {
  // Aggregately construct the arrangement of five line segments.
  Segment segments[] = {Segment(Point(1, 0), Point(2, 4)),
                        Segment(Point(5, 0), Point(5, 5)),
                        Segment(Point(1, 0), Point(5, 3)),
                        Segment(Point(0, 2), Point(6, 0)),
                        Segment(Point(3, 0), Point(5, 5))};
  Arrangement arr;
  insert(arr, segments, segments + sizeof(segments)/sizeof(Segment));
  print_arrangement_size(arr);
  return 0;
}
