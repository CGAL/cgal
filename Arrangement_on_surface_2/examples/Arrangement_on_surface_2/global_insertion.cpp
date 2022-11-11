//! \file examples/Arrangement_on_surface_2/global_insertion.cpp
// Using the global insertion functions (incremental and aggregated).

#include "arr_exact_construction_segments.h"
#include "arr_print.h"

int main() {
  Segment S1[] = {Segment(Point(1, 3), Point(4, 6)),
                  Segment(Point(1, 3), Point(6, 3)),
                  Segment(Point(1, 3), Point(4, 0)),
                  Segment(Point(4, 6), Point(6, 3)),
                  Segment(Point(4, 0), Point(6, 3))};
  Segment s = Segment(Point(0, 3), Point(4, 3));
  Segment S2[] = {Segment(Point(0, 5), Point(6, 6)),
                  Segment(Point(0, 4), Point(6, 5)),
                  Segment(Point(0, 2), Point(6, 1)),
                  Segment(Point(0, 1), Point(6, 0)),
                  Segment(Point(6, 1), Point(6, 5))};

  Arrangement arr;
  insert_non_intersecting_curves(arr, S1, S1 + sizeof(S1)/sizeof(Segment));
  insert(arr, s);                                     // 1 incremental
  insert(arr, S2, S2 + sizeof(S2)/sizeof(Segment));   // 5 aggregate
  print_arrangement_size(arr);
  return 0;
}
