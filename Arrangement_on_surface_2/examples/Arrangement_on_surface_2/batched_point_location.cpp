//! \file examples/Arrangement_on_surface_2/batched_point_location.cpp
// Answering a batched point-location query.

#include <list>

#include <CGAL/basic.h>
#include <CGAL/Arr_batched_point_location.h>

#include "arr_inexact_construction_segments.h"
#include "point_location_utils.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
typedef Traits::Point_2                                 Point;
typedef CGAL::Arrangement_2<Traits>                     Arrangement;
typedef CGAL::Arr_point_location_result<Arrangement>    Point_location_result;
typedef std::pair<Point, Point_location_result::Type>   Query_result;

typedef Arrangement::Vertex_const_handle                Vertex_const_handle;
typedef Arrangement::Halfedge_const_handle              Halfedge_const_handle;
typedef Arrangement::Face_const_handle                  Face_const_handle;

int main()
{
  // Construct the arrangement.
  Arrangement arr;
  construct_segments_arr(arr);

  // Perform a batched point-location query.
  std::list<Point> points;
  points.push_back(Point(1, 4));
  points.push_back(Point(4, 3));
  points.push_back(Point(6, 3));
  points.push_back(Point(3, 2));
  points.push_back(Point(5, 2));
  points.push_back(Point(1, 0));
  std::list<Query_result> results;
  locate(arr, points.begin(), points.end(), std::back_inserter(results));

  // Print the results.
  for (auto it = results.begin(); it != results.end(); ++it)
    print_point_location<Arrangement>(it->first, it->second);
  return 0;
}
