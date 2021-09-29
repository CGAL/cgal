//! \file examples/Arrangement_on_surface_2/dual_lines.cpp
// Checking whether there are three collinear points in a given input set
// using the arrangement of the dual lines.

#include <cstdlib>

#include "arr_linear.h"
#include "read_objects.h"

int main(int argc, char* argv[]) {
  // Get the name of the input file from the command line, or use the default
  // points.dat file if no command-line parameters are given.
  const char* filename = (argc > 1) ? argv[1] : "points.dat";

  // Open the input file.
  std::ifstream in_file(filename);

  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << "!\n";
    return 1;
  }

  // Read the points from the file, and construct their dual lines.
  // The input file format should be (all coordinate values are integers):
  // <n>                                 // number of point.
  // <x_1> <y_1>                         // point #1.
  // <x_2> <y_2>                         // point #2.
  //   :      :       :      :
  // <x_n> <y_n>                         // point #n.
  std::vector<Point> points;
  read_objects<Point>(filename, std::back_inserter(points));
  std::list<X_monotone_curve> dual_lines;
  for (const auto& p : points) dual_lines.push_back(Line(p.x(), -1, -p.y()));

  // Construct the dual arrangement by aggregately inserting the lines.
  Arrangement arr;
  insert(arr, dual_lines.begin(), dual_lines.end());
  std::cout << "The dual arrangement size:\n"
            << "V = " << arr.number_of_vertices()
            << " (+ " << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces()
            << " (" << arr.number_of_unbounded_faces()
            << " unbounded)\n";

  // Look for a vertex whose degree is greater than 4.
  bool found_collinear = false;
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->degree() > 4) {
      found_collinear = true;
      break;
    }
  }
  if (found_collinear)
    std::cout << "Found at least three collinear points in the input set.\n";
  else
    std::cout << "No three collinear points are found in the input set.\n";

  // Pick two points from the input set, compute their midpoint and insert
  // its dual line into the arrangement.
  Kernel ker;
  const auto n = points.size();
  const auto k1 = std::rand() % n, k2 = (k1 + 1) % n;
  Point p_mid = ker.construct_midpoint_2_object()(points[k1], points[k2]);
  X_monotone_curve dual_p_mid = Line(p_mid.x(), -1, -p_mid.y());
  insert(arr, dual_p_mid);

  // Make sure that we now have three collinear points.
  found_collinear = false;
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->degree() > 4) {
      found_collinear = true;
      break;
    }
  }
  CGAL_assertion(found_collinear);

  return (0);
}
