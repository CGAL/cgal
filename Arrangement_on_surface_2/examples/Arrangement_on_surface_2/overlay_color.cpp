//! \file examples/Arrangement_on_surface_2/overlay_color.cpp
// The overlay of two arrangement with extended dcel structures

#include <cassert>

#include <CGAL/basic.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

#include "arr_exact_construction_segments.h"
#include "Overlay_color_traits.h"

using Color = unsigned int;
using Dcel = CGAL::Arr_extended_dcel<Traits, Color, Color, Color>;
using Ex_arrangement = CGAL::Arrangement_2<Traits, Dcel>;

int main() {
  const Color vcol1(0x00000080), hcol1(0x000000ff), fcol1(0x00ccccff);
  const Color vcol2(0x00800000), hcol2(0x00ff0000), fcol2(0x00ffcccc);

  // Construct the first arrangement and assign colors to its features.
  Ex_arrangement arr1;
  insert_non_intersecting_curve(arr1, Segment(Point(0, 0), Point(4, 0)));
  insert_non_intersecting_curve(arr1, Segment(Point(0, 2), Point(4, 2)));
  insert_non_intersecting_curve(arr1, Segment(Point(0, 4), Point(4, 4)));
  insert(arr1, Segment(Point(0, 0), Point(0, 4)));
  insert(arr1, Segment(Point(2, 0), Point(2, 4)));
  insert(arr1, Segment(Point(4, 0), Point(4, 4)));
  assert(arr1.number_of_faces() == 5);
  for (auto vit = arr1.vertices_begin(); vit != arr1.vertices_end(); ++vit)
    vit->set_data(vcol1);
  for (auto hit = arr1.halfedges_begin(); hit != arr1.halfedges_end(); ++hit)
    hit->set_data(hcol1);
  for (auto fit = arr1.faces_begin(); fit != arr1.faces_end(); ++fit)
    fit->set_data(fcol1);

  // Construct the second arrangement and assign colors to its features.
  Ex_arrangement  arr2;
  insert_non_intersecting_curve(arr2, Segment(Point(0, 0), Point(6, 0)));
  insert_non_intersecting_curve(arr2, Segment(Point(0, 3), Point(6, 3)));
  insert_non_intersecting_curve(arr2, Segment(Point(0, 6), Point(6, 6)));
  insert(arr2, Segment(Point(0, 0), Point(0, 6)));
  insert(arr2, Segment(Point(3, 0), Point(3, 6)));
  insert(arr2, Segment(Point(6, 0), Point(6, 6)));
  assert(arr2.number_of_faces() == 5);
  for (auto vit = arr2.vertices_begin(); vit != arr2.vertices_end(); ++vit)
    vit->set_data(vcol2);
  for (auto hit = arr2.halfedges_begin(); hit != arr2.halfedges_end(); ++hit)
    hit->set_data(hcol2);
  for (auto fit = arr2.faces_begin(); fit != arr2.faces_end(); ++fit)
    fit->set_data(fcol2);

  // Compute the overlay of the two arrangements, while blending the colors
  // of their features.
  Ex_arrangement ovl_arr;
  Overlay_color_traits<Ex_arrangement> overlay_traits;
  CGAL::overlay(arr1, arr2, ovl_arr, overlay_traits);

  // Print the overlay-arrangement vertices and their colors.
  for (auto vit = ovl_arr.vertices_begin(); vit != ovl_arr.vertices_end(); ++vit)
    std::cout << vit->point() << ": 0x" << std::hex << std::setfill('0')
              << std::setw(6) << vit->data() << std::endl;
  return 0;
}
