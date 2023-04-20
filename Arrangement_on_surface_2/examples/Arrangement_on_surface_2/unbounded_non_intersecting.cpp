//! \file examples/Arrangement_on_surface_2/unbounded_non_intersecting.cpp
// Constructing an arrangement of unbounded linear objects using the insertion
// function for non-intersecting curves.

#include <cassert>

#include "arr_linear.h"
#include "arr_print.h"

int main() {
  Arrangement arr;

  // Insert a line in the (currently single) unbounded face of the arrangement;
  // then, insert a point that lies on the line splitting it into two.
  X_monotone_curve c1 = Line(Point(-1, 0), Point(1, 0));
  arr.insert_in_face_interior(c1, arr.unbounded_face());
  Vertex_handle v = insert_point(arr, Point(0,0));
  assert(! v->is_at_open_boundary());

  // Add two more rays using the specialized insertion functions.
  arr.insert_from_right_vertex(Ray(Point(0, 0), Point(-1, 1)), v); // c2
  arr.insert_from_left_vertex(Ray(Point(0, 0), Point(1, 1)), v);   // c3

  // Insert three more interior-disjoint rays, c4, c5, and c6.
  insert_non_intersecting_curve(arr, Ray(Point(0, -1), Point(-2, -2)));
  insert_non_intersecting_curve(arr, Ray(Point(0, -1), Point(2, -2)));
  insert_non_intersecting_curve(arr, Ray(Point(0, 0), Point(0, 1)));

  print_unbounded_arrangement_size(arr);

  // Print the outer CCBs of the unbounded faces.
  int k = 1;
  for (auto it = arr.unbounded_faces_begin(); it != arr.unbounded_faces_end();
       ++it)
  {
    std::cout << "Face no. " << k++ << "(" << it->is_unbounded() << ","
              << it->number_of_holes() << ")" << ": ";
    Arrangement::Ccb_halfedge_const_circulator first = it->outer_ccb();
    auto curr = first;
    if (! curr->source()->is_at_open_boundary())
      std::cout << "(" << curr->source()->point() << ")";

    do {
      Arrangement::Halfedge_const_handle e = curr;
      if (! e->is_fictitious()) std::cout << " [" << e->curve() << "] ";
      else std::cout << " [ ... ] ";

      if (! e->target()->is_at_open_boundary())
        std::cout << "(" << e->target()->point() << ")";
    } while (++curr != first);
    std::cout << std::endl;
  }
  return 0;
}
