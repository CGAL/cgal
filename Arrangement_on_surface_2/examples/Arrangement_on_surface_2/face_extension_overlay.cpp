//! \file examples/Arrangement_on_surface_2/face_extension_overlay.cpp
// A face overlay of two arrangements with extended face records.

#include <cassert>

#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

#include "arr_exact_construction_segments.h"

using Dcel = CGAL::Arr_face_extended_dcel<Traits, bool>;
using Ex_arrangement = CGAL::Arrangement_2<Traits, Dcel>;
using Overlay_traits =
  CGAL::Arr_face_overlay_traits<Ex_arrangement, Ex_arrangement,
                                Ex_arrangement, std::logical_and<bool>>;

int main() {
  // Construct the first arrangement, containing a square-shaped face.
  Ex_arrangement arr1;
  insert_non_intersecting_curve(arr1, Segment(Point(2, 2), Point(6, 2)));
  insert_non_intersecting_curve(arr1, Segment(Point(6, 2), Point(6, 6)));
  insert_non_intersecting_curve(arr1, Segment(Point(6, 6), Point(2, 6)));
  insert_non_intersecting_curve(arr1, Segment(Point(2, 6), Point(2, 2)));
  // 2 because the bounded and the unbounded one
  assert(arr1.number_of_faces() == 2);

  // Mark just the bounded face.
  for (auto fit = arr1.faces_begin(); fit != arr1.faces_end(); ++fit)
    fit->set_data(fit != arr1.unbounded_face());

  // Construct the second arrangement, containing a rhombus-shaped face.
  Ex_arrangement arr2;
  insert_non_intersecting_curve(arr2, Segment(Point(4, 1), Point(7, 4)));
  insert_non_intersecting_curve(arr2, Segment(Point(7, 4), Point(4, 7)));
  insert_non_intersecting_curve(arr2, Segment(Point(4, 7), Point(1, 4)));
  insert_non_intersecting_curve(arr2, Segment(Point(1, 4), Point(4, 1)));

  for (auto fit = arr2.faces_begin(); fit != arr2.faces_end(); ++fit)
    fit->set_data(fit != arr2.unbounded_face());    // mark the bounded face.

  // Compute the overlay of the two arrangements, marking only the faces that
  // are intersections of two marked faces in arr1 and arr2, respectively.
  Ex_arrangement overlay_arr;
  Overlay_traits overlay_traits;
  CGAL::overlay(arr1, arr2, overlay_arr, overlay_traits);

  // Go over the faces of the resulting arrangement and print the marked ones.
  std::cout << "The intersection is: ";
  for (auto fit = overlay_arr.faces_begin(); fit != overlay_arr.faces_end();
       ++fit)
  {
    if (! fit->data()) continue;
    Ex_arrangement::Ccb_halfedge_circulator curr = fit->outer_ccb();
    std::cout << curr->source()->point();
    do std::cout << " --> " << curr->target()->point();
    while (++curr != fit->outer_ccb());
    std::cout << std::endl;
  }
  return 0;
}
