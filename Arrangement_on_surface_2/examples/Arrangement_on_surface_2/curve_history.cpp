//! \file examples/Arrangement_on_surface_2/curve_history.cpp
// Constructing an arrangement with curve history.

#include <CGAL/basic.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>

#include "arr_exact_construction_segments.h"
#include "point_location_utils.h"

typedef CGAL::Arrangement_with_history_2<Traits>              Arr_with_hist;
typedef Arr_with_hist::Curve_handle                           Curve_handle;
typedef CGAL::Arr_trapezoid_ric_point_location<Arr_with_hist> Point_location;

int main() {
  // Insert 3 curves incrementally.
  Arr_with_hist arr;
  insert(arr, Segment(Point(0, 3), Point(4, 3)));
  insert(arr, Segment(Point(3, 2), Point(3, 5)));
  insert(arr, Segment(Point(2, 3), Point(5, 3)));

  // Insert three additional segments aggregately.
  Segment segs[] = {Segment(Point(2, 6), Point(7, 1)),
                    Segment(Point(0, 0), Point(2, 6)),
                    Segment(Point(3, 4), Point(6, 4))};
  insert(arr, segs, segs + sizeof(segs)/sizeof(Segment));

  // Print out the curves and the number of edges each one induces.
  std::cout << "The arrangement contains "
            << arr.number_of_curves() << " curves:\n";
  for (auto cit = arr.curves_begin(); cit != arr.curves_end(); ++cit)
    std::cout << "Curve [" << *cit << "] induces "
              << arr.number_of_induced_edges(cit) << " edges.\n";

  // Print the arrangement edges along with the list of curves that
  // induce each edge.
  std::cout << "The arrangement comprises "
            << arr.number_of_edges() << " edges:\n";
  for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    std::cout << "[" << eit->curve() << "]. Originating curves: ";
    for (auto ocit = arr.originating_curves_begin(eit);
         ocit != arr.originating_curves_end(eit); ++ocit)
      std::cout << " [" << *ocit << "]" << std::flush;
    std::cout << std::endl;
  }

  // Perform some point-location queries.
  Point_location pl(arr);
  locate_point(pl, Point(4, 6));      // q1
  locate_point(pl, Point(6, 2));      // q2
  locate_point(pl, Point(2, 4));      // q3

  return 0;
}
