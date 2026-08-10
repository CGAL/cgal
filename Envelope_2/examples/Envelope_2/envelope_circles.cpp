//! \file examples/Example_2/envelope_circles.cpp
// Constructing the envelopes of a set of circles using the circle-segment
// traits.

#include <cassert>

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Envelope_2/Envelope_diagram_1.h>
#include <CGAL/Envelope_2/envelope_2.h>

using Number_type = CGAL::Exact_rational;
using Kernel = CGAL::Cartesian<Number_type>;
using Kernel_point_2 = Kernel::Point_2;
using Circle_2 = Kernel::Circle_2;
using Traits_2 = CGAL::Arr_circle_segment_traits_2<Kernel>;
using Curve_2 = Traits_2::Curve_2;
using Diagram_1 = CGAL::Envelope_2::Envelope_diagram_1<Traits_2>;

/*! Print the given envelope diagram. */
void print_diagram (const Diagram_1& diag) {
  Diagram_1::Edge_const_descriptor e = diag.leftmost();
  while (e != diag.rightmost()) {
    std::cout << "Edge: ";
    if (! diag.empty_edge_curves(e)) {
      Circle_2 circ = diag.edge_curve(e).supporting_circle();
      std::cout << " (x - " << CGAL::to_double(circ.center().x()) << ")^2 +"
                << " (y - " << CGAL::to_double(circ.center().y()) << ")^2 = "
                << CGAL::to_double(circ.squared_radius()) << std::endl;
    }
    else std::cout << " [empty]" << std::endl;
    auto v = diag.right_vertex(e);
    std::cout << "Vertex (" << CGAL::to_double(diag.point(v).x()) << ' '
              << CGAL::to_double(diag.point(v).y()) << ')' << std::endl;
    e = diag.right_edge(v);
  }
  assert(diag.empty_edge_curves(e));
  std::cout << "Edge: [empty]" << std::endl;
}

/*! The main function. */
int main() {
  // Create four input circles.
  Curve_2 circles[4];

  circles[0] = Circle_2(Kernel_point_2(1, 3), CGAL::square(2));
  circles[1] = Circle_2(Kernel_point_2(4, 5), CGAL::square(4));
  circles[2] = Circle_2(Kernel_point_2(5, 1), CGAL::square(1));
  circles[3] = Circle_2(Kernel_point_2(6, 7), CGAL::square(2));

  // Compute the minimization diagram that represents their lower envelope.
  Diagram_1 min_diag;

  CGAL::Envelope_2::lower_envelope_2(&(circles[0]), &(circles[4]), min_diag);
  print_diagram(min_diag);

  // Compute the maximization diagram that represents the upper envelope.
  Diagram_1 max_diag;

  CGAL::Envelope_2::upper_envelope_2(&(circles[0]), &(circles[4]), max_diag);
  print_diagram(max_diag);

  return 0;
}
