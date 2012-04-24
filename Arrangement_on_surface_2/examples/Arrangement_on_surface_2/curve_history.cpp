//! \file examples/Arrangement_on_surface_2/curve_history.cpp
// Constructing an arrangement with curve history.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>

#include <CGAL/Arrangement_on_surface_with_history_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_simple_point_location.h>

#include "point_location_utils.h"

typedef CGAL::Cartesian<Number_type>                      Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                Traits_2;
typedef Traits_2::Point_2                                 Point_2;
typedef Traits_2::Curve_2                                 Segment_2;
typedef CGAL::Arrangement_with_history_2<Traits_2>        Arr_with_hist_2;
typedef Arr_with_hist_2::Curve_handle                     Curve_handle;
typedef CGAL::Arr_simple_point_location<Arr_with_hist_2>  Point_location;

int main()
{
  Arr_with_hist_2   arr;

  // Insert s1, s2 and s3 incrementally:
  Segment_2 s1(Point_2(0, 3), Point_2(4, 3));
  insert(arr, s1);
  Segment_2 s2(Point_2(3, 2), Point_2(3, 5));
  insert(arr, s2);
  Segment_2 s3(Point_2(2, 3), Point_2(5, 3));
  insert(arr, s3);

  // Insert three additional segments aggregately:
  Segment_2 segs[3];
  segs[0] = Segment_2(Point_2(2, 6), Point_2(7, 1));
  segs[1] = Segment_2(Point_2(0, 0), Point_2(2, 6));
  segs[2] = Segment_2(Point_2(3, 4), Point_2(6, 4));
  insert(arr, segs, segs + 3);

  // Print out the curves and the number of edges each one induces.
  Arr_with_hist_2::Curve_iterator            cit;
  std::cout << "The arrangement contains "
            << arr.number_of_curves() << " curves:" << std::endl;
  for (cit = arr.curves_begin(); cit != arr.curves_end(); ++cit)
    std::cout << "Curve [" << *cit << "] induces "
              << arr.number_of_induced_edges(cit) << " edges." << std::endl; 

  // Print the arrangement edges, along with the list of curves that
  // induce each edge.
  Arr_with_hist_2::Edge_iterator                  eit;
  Arr_with_hist_2::Originating_curve_iterator     ocit;
  std::cout << "The arrangement is comprised of "
            << arr.number_of_edges() << " edges:" << std::endl;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    std::cout << "[" << eit->curve() << "]. Originating curves: ";
    for (ocit = arr.originating_curves_begin(eit);
         ocit != arr.originating_curves_end(eit); ++ocit)
      std::cout << " [" << *ocit << "]" << std::flush;
    std::cout << std::endl;
  }

  // Perform some point-location queries:
  Point_location   pl(arr);

  Point_2          p1(4, 6);
  point_location_query(pl, p1);
  Point_2          p2(6, 2);
  point_location_query(pl, p2);
  Point_2          p3(2, 4);
  point_location_query(pl, p3);

  return 0;
}
