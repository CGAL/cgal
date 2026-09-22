//! \file examples/Arrangement_on_surface_2/edge_manipulation_curve_history.cpp
// Removing curves and manipulating edges in an arrangement with history.

#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

#include "arr_circular.h"
#include "arr_print.h"

using Arr_with_hist = CGAL::Arrangement_with_history_2<Traits>;
using Curve_handle = Arr_with_hist::Curve_handle;
using Point_location = CGAL::Arr_walk_along_line_point_location<Arr_with_hist>;

int main() {
  // Construct an arrangement containing nine circles: C[0] of radius 2 and
  // C[1], ..., C[8] of radius 1.
  const Number_type _7_halves = Number_type(7) / Number_type(2);
  Curve C[9];

  C[0] = Circle(Kernel::Point_2(_7_halves, _7_halves), 4, CGAL::CLOCKWISE);
  C[1] = Circle(Kernel::Point_2(_7_halves, 6), 1, CGAL::CLOCKWISE);
  C[2] = Circle(Kernel::Point_2(5, 6), 1, CGAL::CLOCKWISE);
  C[3] = Circle(Kernel::Point_2(6, _7_halves), 1, CGAL::CLOCKWISE);
  C[4] = Circle(Kernel::Point_2(5, 2), 1, CGAL::CLOCKWISE);
  C[5] = Circle(Kernel::Point_2(_7_halves, 1), 1, CGAL::CLOCKWISE);
  C[6] = Circle(Kernel::Point_2(2, 2), 1, CGAL::CLOCKWISE);
  C[7] = Circle(Kernel::Point_2(1, _7_halves), 1, CGAL::CLOCKWISE);
  C[8] = Circle(Kernel::Point_2(2, 5), 1, CGAL::CLOCKWISE);

  Arr_with_hist arr;
  Curve_handle handles[9];
  for (size_t k = 0; k < 9; ++k) handles[k] = insert(arr, C[k]);

  std::cout << "The initial arrangement size:\n";
  print_arrangement_size(arr);

  // Remove the large circle C[0].
  std::cout << "Removing C[0]: ";
  std::cout << remove_curve(arr, handles[0])
            << " edges have been removed.\n";
  print_arrangement_size(arr);

  // Locate the point q, which should be on an edge e.
  Point_location pl(arr);
  const Point q{_7_halves, 7};
  Point_location::result_type obj = pl.locate(q);
  auto* e = std::get_if<Arr_with_hist::Halfedge_const_handle>(&obj);

  // Split the edge e to two edges e1 and e2;
  auto e1 = arr.split_edge(arr.non_const_handle(*e), q);
  auto e2 = e1->next();
  std::cout << "After edge split:\n";
  print_arrangement_size(arr);

  // Merge back the two split edges.
  arr.merge_edge(e1, e2);
  std::cout << "After edge merge:\n";
  print_arrangement_size(arr);
  return 0;
}
