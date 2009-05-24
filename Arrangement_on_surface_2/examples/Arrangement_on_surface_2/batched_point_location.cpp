//! \file examples/Arrangement_on_surface_2/batched_point_location.cpp
// Answering a batched point-location query.

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_batched_point_location.h>
#include <list>

#include "point_location_utils.h"

typedef CGAL::MP_Float                                Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef std::pair<Point_2, CGAL::Object>              Query_result;

int main ()
{
  // Construct the arrangement.
  Arrangement_2    arr;

  construct_segments_arr (arr);

  // Perform a batched point-location query.
  std::list<Point_2>       query_points;
  std::list<Query_result>  results;

  query_points.push_back (Point_2 (1, 4));
  query_points.push_back (Point_2 (4, 3));
  query_points.push_back (Point_2 (6, 3));
  query_points.push_back (Point_2 (3, 2));
  query_points.push_back (Point_2 (5, 2));
  query_points.push_back (Point_2 (1, 0));

  locate (arr, query_points.begin(), query_points.end(),
	  std::back_inserter (results));

  // Print the results.
  std::list<Query_result>::const_iterator  res_iter;
  Arrangement_2::Vertex_const_handle       v;
  Arrangement_2::Halfedge_const_handle     e;
  Arrangement_2::Face_const_handle         f;

  for (res_iter = results.begin(); res_iter != results.end(); ++res_iter) {
    std::cout << "The point (" << res_iter->first << ") is located ";
    if (CGAL::assign (f, res_iter->second)) {
      // The current qeury point is located inside a face:
      if (f->is_unbounded())
        std::cout << "inside the unbounded face." << std::endl;
      else
        std::cout << "inside a bounded face." << std::endl;
    }
    else if (CGAL::assign (e, res_iter->second)) {
      // The current qeury point is located on an edge:
      std::cout << "on an edge: " << e->curve() << std::endl;
    }
    else if (CGAL::assign (v, res_iter->second)) {
      // The current qeury point is located on a vertex:
      if (v->is_isolated())
        std::cout << "on an isolated vertex: " << v->point() << std::endl;
      else
        std::cout << "on a vertex: " << v->point() << std::endl;
    }
  }

  return 0;
}
