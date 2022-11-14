//! \file examples/Arrangement_2/io_unbounded.cpp
// Using the I/O operators with an arrangement of unbounded curves.

#include <list>
#include <fstream>

#include <CGAL/basic.h>
#include <CGAL/IO/Arr_iostream.h>

#include "arr_linear.h"

int main() {
  // Construct an arrangement of five linear objects.
  Arrangement arr;
  std::list<X_monotone_curve> curves;

  curves.push_back(Line(Point(0, 0), Point(2, 1)));
  curves.push_back(Line(Point(0, 0), Point(2, -1)));
  curves.push_back(Line(Point(-1, 0), Point(-1, 1)));
  curves.push_back(Ray(Point(2, 3), Point(2, 4)));
  curves.push_back(Segment(Point(0, 1), Point(0, 2)));

  insert(arr, curves.begin(), curves.end());

  // Print out the size of the resulting arrangement.
  std::cout << "Writing an arrangement of size:\n"
            << "   V = " << arr.number_of_vertices()
            << " (plus " << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces()
            << " (" << arr.number_of_unbounded_faces() << " unbounded)\n\n";

  // Write the arrangement to a file.
  std::ofstream out_file("arr_ex_io_unbounded.dat");

  out_file << arr;
  out_file.close();

  // Read the arrangement from the file.
  Arrangement arr2;
  std::ifstream in_file("arr_ex_io_unbounded.dat");

  in_file >> arr2;
  in_file.close();

  std::cout << "Read an arrangement of size:\n"
            << "   V = " << arr2.number_of_vertices()
            << " (plus " << arr2.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr2.number_of_edges()
            << ",  F = " << arr2.number_of_faces()
            << " (" << arr2.number_of_unbounded_faces() << " unbounded)\n\n";

  return 0;
}
