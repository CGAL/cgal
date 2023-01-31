//! \file examples/Arrangement_on_surface_2/dcel_extension.cpp
// Extending all DCEL records (vertices, edges and faces).

#include <CGAL/basic.h>
#include <CGAL/Arr_extended_dcel.h>

#include "arr_exact_construction_segments.h"

enum Color {BLUE, RED, WHITE};

typedef CGAL::Arr_extended_dcel<Traits, Color, bool, size_t> Dcel;
typedef CGAL::Arrangement_2<Traits, Dcel>                    Ex_arrangement;

int main() {
  // Construct the arrangement containing two intersecting triangles.
  Traits traits;
  Ex_arrangement arr(&traits);
  insert_non_intersecting_curve(arr, Segment(Point(4, 1), Point(7, 6)));
  insert_non_intersecting_curve(arr, Segment(Point(1, 6), Point(7, 6)));
  insert_non_intersecting_curve(arr, Segment(Point(4, 1), Point(1, 6)));
  insert(arr, Segment(Point(1, 3), Point(7, 3)));
  insert(arr, Segment(Point(1, 3), Point(4, 8)));
  insert(arr, Segment(Point(4, 8), Point(7, 3)));
  insert_point(arr, Point(4, 4.5));

  // Go over all arrangement edges and set their flags.
  // Recall that the value type of the edge iterator is the halfedge type.
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    auto degree = vit->degree();
    vit->set_data((degree == 0) ? BLUE : ((degree <= 2) ? RED : WHITE));
  }

  auto equal = traits.equal_2_object();
  for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    // Check whether the halfedge has the same direction as its segment.
    bool flag = equal(eit->source()->point(),eit->curve().source());
    eit->set_data(flag);
    eit->twin()->set_data(!flag);
  }

  // Store the size of the outer boundary of every face of the arrangement.
  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    size_t boundary_size = 0;
    if (! fit->is_unbounded()) {
      Ex_arrangement::Ccb_halfedge_circulator curr = fit->outer_ccb();
      boundary_size = std::distance(++curr, fit->outer_ccb())+1;
    }
    fit->set_data(boundary_size);
  }

  // Copy the arrangement and print the vertices along with their colors.
  Ex_arrangement arr2 = arr;

  std::cout << "The arrangement vertices:\n";
  for (auto vit = arr2.vertices_begin(); vit != arr2.vertices_end(); ++vit) {
    std::cout << '(' << vit->point() << ") - ";
    switch (vit->data()) {
      case BLUE  : std::cout << "BLUE.\n"; break;
      case RED   : std::cout << "RED.\n"; break;
      case WHITE : std::cout << "WHITE.\n"; break;
    }
  }

  std::cout << "The arrangement outer-boundary sizes:";
  for (auto fit = arr2.faces_begin(); fit != arr2.faces_end(); ++fit)
    std::cout << " " << fit->data();
  std::cout << std::endl;
  return 0;
}
