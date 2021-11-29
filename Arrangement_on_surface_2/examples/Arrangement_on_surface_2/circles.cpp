//! \file examples/Arrangement_on_surface_2/circles.cpp
// Constructing an arrangement of circles using the circle-segment traits.

#include "arr_circular.h"

int main() {
  Arrangement arr;

  // Create a circle centered at the origin with radius 5 (C1).
  insert(arr, Curve(Circle(Rational_point(0, 0), Number_type(25))));

  // Create a circle centered at (7,7) with radius 5 (C2).
  insert(arr, Curve(Circle(Rational_point(7, 7), Number_type(25))));

  // Create a circle centered at (4,-0.5) with radius 3.5 (= 7/2) (C3).
  Rational_point c3(4, Number_type(-1) / Number_type(2));
  insert(arr, Curve(Circle(c3, Number_type(49) / Number_type(4))));

  // Locate the vertex with maximal degree.
  auto vit = arr.vertices_begin();
  Arrangement::Vertex_const_handle v_max = vit;
  for (++vit; vit != arr.vertices_end(); ++vit)
    if (vit->degree() > v_max->degree()) v_max = vit;

  // Locate the vertex with maximum degree.

  std::cout << "The vertex with maximal degree in the arrangement is: "
            << "v_max = (" << v_max->point() << ") "
            << "with degree " << v_max->degree() << "." << std::endl;
  return 0;
}
