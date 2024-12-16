//! \file examples/Arrangement_on_surface_2/batched_point_location.cpp
// Answering a batched point-location query.

#include <list>
#include <vector>

#include <CGAL/basic.h>
#include <CGAL/Arr_batched_point_location.h>

#include "arr_inexact_construction_segments.h"
#include "point_location_utils.h"

using Point_location_result = CGAL::Arr_point_location_result<Arrangement>;
using Query_result = std::pair<Point, Point_location_result::Type>;

int main() {
  // Construct the arrangement.
  Arrangement arr;
  construct_segments_arr(arr);

  // Perform a batched point-location query.
  std::vector<Point> points = {
    Point(1, 4), Point(4, 3), Point(6, 3), Point(3, 2), Point(5, 2), Point(1, 0)
  };
  std::list<Query_result> results;
  CGAL::locate(arr, points.begin(), points.end(), std::back_inserter(results));

  // Print the results.
  for (auto it = results.begin(); it != results.end(); ++it)
    print_point_location<Arrangement>(it->first, it->second);
  return 0;
}
