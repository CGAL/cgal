//! \file examples/Example_2/envelope_segments.cpp
// Constructing the lower envelope of a set of segments using the new Envelope_diagram_1 API.

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Envelope_diagram_1.h>
#include <CGAL/envelope_2.h>

#include <list>
#include <iostream>
#include <cassert>

using Number_type = CGAL::Exact_rational;
using Kernel = CGAL::Cartesian<Number_type>;
using Segment_traits_2 = CGAL::Arr_segment_traits_2<Kernel>;
using Segment_2 = Segment_traits_2::X_monotone_curve_2;
using Traits_2 = CGAL::Arr_curve_data_traits_2<Segment_traits_2, char>;
using Point_2 = Traits_2::Point_2;
using Labeled_segment_2 = Traits_2::X_monotone_curve_2;
using Diagram_1 = CGAL::Envelope_diagram_1<Traits_2>;

int main() {
  // Construct the input segments and label them 'A' ... 'H'.
  std::list<Labeled_segment_2> segments;

  segments.push_back(Labeled_segment_2(Segment_2(Point_2(0, 1), Point_2(2, 3)), 'A'));
  segments.push_back(Labeled_segment_2(Segment_2(Point_2(1, 2), Point_2(4, 5)), 'B'));
  segments.push_back(Labeled_segment_2(Segment_2(Point_2(1, 5), Point_2(7, 2)), 'C'));
  segments.push_back(Labeled_segment_2(Segment_2(Point_2(4, 2), Point_2(6, 4)), 'D'));
  segments.push_back(Labeled_segment_2(Segment_2(Point_2(8, 3), Point_2(8, 6)), 'E'));
  segments.push_back(Labeled_segment_2(Segment_2(Point_2(9, 2), Point_2(12, 4)), 'F'));
  segments.push_back(Labeled_segment_2(Segment_2(Point_2(10, 2), Point_2(12, 1)), 'G'));
  segments.push_back(Labeled_segment_2(Segment_2(Point_2(11, 0), Point_2(11, 5)), 'H'));

  // Compute the minimization diagram that represents their lower envelope.
  Diagram_1 min_diag;
  CGAL::lower_envelope_x_monotone_2(segments.begin(), segments.end(), min_diag);

  // Print the minimization diagram using the new API functions on min_diag.
  Diagram_1::Edge_const_descriptor e = min_diag.leftmost();
  while (e != min_diag.rightmost()) {
    std::cout << "Edge:";
    if (! min_diag.empty_edge_curves(e)) {
      for (const auto& cv : min_diag.edge_curves(e)) std::cout << ' ' << cv.data();
      std::cout << std::endl;
    }
    else std::cout << " [empty]" << std::endl;
    auto v = min_diag.right_vertex(e);
    std::cout << "Vertex (" << min_diag.point(v) << "):";
    for (const auto& cv : min_diag.vertex_curves(v)) std::cout << ' ' << cv.data();
    std::cout << std::endl;
    e = min_diag.right_edge(v);
  }
  assert(min_diag.empty_edge_curves(e));
  std::cout << "Edge: [empty]" << std::endl;
  return 0;
}
