//! \file examples/Arrangement_on_surface_2/io_curve_history.cpp
// Using the arrangement-with-history I/O operators.

#include <fstream>

#include <CGAL/basic.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/IO/Arr_with_history_iostream.h>

#include "arr_exact_construction_segments.h"
#include "arr_print.h"

using Arr_with_hist = CGAL::Arrangement_with_history_2<Traits>;

int main() {
  // Insert six additional segments aggregately:
  Segment segs[6];
  segs[0] = Segment(Point(2, 6), Point(7, 1));
  segs[1] = Segment(Point(3, 2), Point(3, 5));
  segs[2] = Segment(Point(2, 3), Point(5, 3));
  segs[3] = Segment(Point(2, 6), Point(7, 1));
  segs[4] = Segment(Point(0, 0), Point(2, 6));
  segs[5] = Segment(Point(3, 4), Point(6, 4));

  Arr_with_hist arr1;
  insert(arr1, segs, segs + 6);
  std::cout << "Writing an arrangement of "
            << arr1.number_of_curves() << " input segments:\n";
  print_arrangement_size(arr1);

  // Write the arrangement to a file.
  std::ofstream out_file("arr_ex_io_hist.dat");
  out_file << arr1;
  out_file.close();

  // Read the arrangement from the file.
  Arr_with_hist arr2;
  std::ifstream in_file("arr_ex_io_hist.dat");
  in_file >> arr2;
  in_file.close();
  std::cout << "Read an arrangement of "
            << arr2.number_of_curves() << " input segments:\n";
  print_arrangement_size(arr2);

  return 0;
}
