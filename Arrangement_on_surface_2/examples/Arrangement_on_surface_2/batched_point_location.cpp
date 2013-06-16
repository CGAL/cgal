//! \file examples/Arrangement_on_surface_2/batched_point_location.cpp
// Answering a batched point-location query.

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_batched_point_location.h>
#include <list>

#include "point_location_utils.h"

typedef CGAL::MP_Float                                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;
typedef CGAL::Arr_point_location_result<Arrangement_2>  Point_location_result;
typedef std::pair<Point_2, Point_location_result::Type> Query_result;

typedef Arrangement_2::Vertex_const_handle              Vertex_const_handle;
typedef Arrangement_2::Halfedge_const_handle            Halfedge_const_handle;
typedef Arrangement_2::Face_const_handle                Face_const_handle;

int main ()
{
  // Construct the arrangement.
  Arrangement_2    arr;
  construct_segments_arr(arr);

  // Perform a batched point-location query.
  std::list<Point_2> points;
  points.push_back(Point_2(1, 4));
  points.push_back(Point_2(4, 3));
  points.push_back(Point_2(6, 3));
  points.push_back(Point_2(3, 2));
  points.push_back(Point_2(5, 2));
  points.push_back(Point_2(1, 0));
  std::list<Query_result>  results;
  locate(arr, points.begin(), points.end(), std::back_inserter(results));

  // Print the results.
  std::list<Query_result>::const_iterator it;
  for (it = results.begin(); it != results.end(); ++it) {
    std::cout << "The point (" << it->first << ") is located ";
    if (const Face_const_handle* f = boost::get<Face_const_handle>(&(it->second)))       // inside a face
      std::cout << "inside "
                << (((*f)->is_unbounded()) ? "the unbounded" : "a bounded")
                << " face." << std::endl;
    else if (const Halfedge_const_handle* e = boost::get<Halfedge_const_handle>(&(it->second))) // on an edge
      std::cout << "on an edge: " << (*e)->curve() << std::endl;
    else if (const Vertex_const_handle*   v = boost::get<Vertex_const_handle>(&(it->second)))  // on a vertex
      std::cout << "on "
                << (((*v)->is_isolated()) ? "an isolated" : "a")
                << " vertex: " << (*v)->point() << std::endl;
  }

  return 0;
}
