//! \file examples/Arrangement_on_surface_2/io.cpp
// Using the arrangement I/O operators.

#include <fstream>

#include <CGAL/IO/Arr_iostream.h>

#include "arr_inexact_construction_segments.h"
#include "arr_print.h"
#include "point_location_utils.h"

int main() {
  // Construct the arrangement.
  Arrangement arr1;
  construct_segments_arr(arr1);
  std::cout << "Writing\n";
  print_arrangement_size(arr1);

  // Write the arrangement to a file.
  std::ofstream out_file("arr_ex_io.dat");
  out_file << arr1;
  out_file.close();

  // Read the arrangement from the file.
  Arrangement arr2;
  std::ifstream in_file("arr_ex_io.dat");
  in_file >> arr2;
  in_file.close();
  std::cout << "Reading\n";
  print_arrangement_size(arr2);

  return 0;
}
