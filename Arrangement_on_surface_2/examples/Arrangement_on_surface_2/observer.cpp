//! \file examples/Arrangement_on_surface_2/observer.cpp
// Using a simple arrangement observer.

#include <CGAL/basic.h>
#include <CGAL/Arr_observer.h>

#include "arr_exact_construction_segments.h"
#include "arr_print.h"

// An observer that receives notifications of face splits and face mergers.
class My_observer : public CGAL::Arr_observer<Arrangement> {
public:
  My_observer(Arrangement& arr) : CGAL::Arr_observer<Arrangement>(arr) {}

  virtual void before_split_face(Face_handle, Halfedge_handle e) {
    std::cout << "-> The insertion of :  [ " << e->curve()
              << " ]  causes a face to split.\n";
  }

  virtual void before_merge_face(Face_handle, Face_handle, Halfedge_handle e) {
    std::cout << "-> The removal of :  [ " << e->curve()
              << " ]  causes two faces to merge.\n";
  }
};

int main() {
  // Construct the arrangement containing one diamond-shaped face.
  Arrangement arr;
  My_observer obs(arr);
  insert_non_intersecting_curve(arr, Segment(Point(-1, 0), Point(0, 1)));
  insert_non_intersecting_curve(arr, Segment(Point(0, 1), Point(1, 0)));
  insert_non_intersecting_curve(arr, Segment(Point(1, 0), Point(0, -1)));
  insert_non_intersecting_curve(arr, Segment(Point(0, -1), Point(-1, 0)));

  // Insert a vertical segment dividing the diamond into two, and a
  // a horizontal segment further dividing the diamond into four.
  Segment s_v(Point(0, -1), Point(0, 1));
  Halfedge_handle e_v = insert_non_intersecting_curve(arr, s_v);
  insert(arr, Segment(Point(-1, 0), Point(1, 0))); /* \label{lst:observer:insertion} */
  print_arrangement_size(arr);

  // Now remove a portion of the vertical segment.
  remove_edge(arr, e_v);            // the observer will make a notification
  print_arrangement_size(arr);

  return 0;
}
