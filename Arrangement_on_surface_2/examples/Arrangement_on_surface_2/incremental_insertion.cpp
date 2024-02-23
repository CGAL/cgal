//! \file examples/Arrangement_on_surface_2/incremental_insertion.cpp
// Using the global incremental insertion functions.

#include <CGAL/basic.h>
#include <CGAL/Arr_naive_point_location.h>

#include "arr_exact_construction_segments.h"
#include "arr_print.h"

using Naive_pl = CGAL::Arr_naive_point_location<Arrangement>;
using Pl_result_type = CGAL::Arr_point_location_result<Arrangement>::Type;

int main() {
  // Construct the arrangement of five intersecting segments.
  Arrangement arr;
  Naive_pl pl(arr);

  Segment s1(Point(1, 0), Point(2, 4));
  Segment s2(Point(5, 0), Point(5, 5));
  Segment s3(Point(1, 0), Point(5, 3));
  Segment s4(Point(0, 2), Point(6, 0));
  Segment s5(Point(3, 0), Point(5, 5));

  auto e = insert_non_intersecting_curve(arr, s1, pl);
  insert_non_intersecting_curve(arr, s2, pl);
  insert(arr, s3, Pl_result_type(e->source()));
  insert(arr, s4, pl);
  insert(arr, s5, pl);
  print_arrangement_size(arr);

  // Perform a point-location query on the resulting arrangement and print
  // the boundary of the face that contains it.
  Point q(4, 1);
  auto obj = pl.locate(q);
  auto* f = std::get_if<Arrangement::Face_const_handle>(&obj);

  std::cout << "The query point (" << q << ") is located in: ";
  print_face<Arrangement>(*f);

  return 0;
}
