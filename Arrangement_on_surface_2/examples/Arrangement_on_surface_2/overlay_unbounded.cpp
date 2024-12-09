//! \file examples/Arrangement_on_surface_2/overlay_unbounded.cpp
// A face overlay of two arrangements with unbounded faces.

#include <string>

#include <CGAL/basic.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

#include "arr_linear.h"

// Define a functor for creating a label from a character and an integer.
struct Overlay_label {
  std::string operator()(char c, unsigned int i) const
  { return c + std::to_string(i); }
};

using Dcel_dlue = CGAL::Arr_face_extended_dcel<Traits, char>;
using Arrangement_blue = CGAL::Arrangement_2<Traits, Dcel_dlue>;
using Dcel_red = CGAL::Arr_face_extended_dcel<Traits, unsigned int>;
using Arrangement_red = CGAL::Arrangement_2<Traits, Dcel_red>;
using Dcel_res = CGAL::Arr_face_extended_dcel<Traits, std::string>;
using Arrangement_res = CGAL::Arrangement_2<Traits, Dcel_res>;
using Overlay_traits =
  CGAL::Arr_face_overlay_traits<Arrangement_blue, Arrangement_red,
                                Arrangement_res, Overlay_label>;

int main() {
  // Construct the first arrangement, induced by two lines y = x and y = -x.
  Arrangement_blue arr1;
  insert(arr1, Line(Point(0, 0), Point(1, 1)));
  insert(arr1, Line(Point(0, 0), Point(1, -1)));

  // Label the four (unbounded) faces of the arrangement as 'A' to 'D' by
  // traversing the faces incident to the halfedges around the single
  // arrangement vertex (0, 0).
  char clabel = 'A';
  auto first = arr1.vertices_begin()->incident_halfedges();
  auto curr = first;
  do curr->face()->set_data(clabel++);
  while (++curr != first);

  // Construct the second arrangement, containing a single square-shaped face.
  Arrangement_red arr2;
  insert_non_intersecting_curve(arr2, Segment(Point(-3, -3), Point(3, -3)));
  insert_non_intersecting_curve(arr2, Segment(Point(3, -3), Point(3, 3)));
  insert_non_intersecting_curve(arr2, Segment(Point(3, 3), Point(-3, 3)));
  insert_non_intersecting_curve(arr2, Segment(Point(-3, 3), Point(-3, -3)));

  // Give the unbounded face the index 1, and the bounded face the index 2.
  for (auto fit = arr2.faces_begin(); fit != arr2.faces_end(); ++fit)
    fit->set_data((fit == arr2.unbounded_face()) ? 1 : 2);

  // Compute the overlay of the two arrangements.
  Arrangement_res overlay_arr;
  Overlay_traits overlay_traits;
  CGAL::overlay(arr1, arr2, overlay_arr, overlay_traits);

  // Go over the faces of the overlay arrangement and print their labels.
  std::cout << "The overlay faces are:\n";
  for (auto res_fit = overlay_arr.faces_begin();
       res_fit != overlay_arr.faces_end(); ++res_fit)
  {
    std::cout << "  " << res_fit->data().c_str() << " ("
              << (res_fit->is_unbounded() ? "unbounded" : "bounded") << ").\n";
  }
  return 0;
}
