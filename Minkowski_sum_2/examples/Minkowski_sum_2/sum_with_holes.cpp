//! \file examples/Minkowski_sum_2/sum_with_holes.cpp
// Computing the Minkowski sum of two non-convex polygons read from a file.

#include <fstream>

#include <CGAL/basic.h>
#include <CGAL/minkowski_sum_2.h>

#include "bops_linear.h"
#include "pgn_print.h"

int main(int argc, char* argv[])
{
  // Open the input file and read the two polygons from it.
  const char* filename = (argc > 1) ? argv[1] : "rooms_star.dat";
  std::ifstream    in_file(filename);
  if (! in_file.is_open()) {
    std::cerr << "Failed to open the input file." << std::endl;
    return -1;
  }
  Polygon_2   P, Q;
  in_file >> P >> Q;
  in_file.close();

  // Compute and print the Minkowski sum.
  Polygon_with_holes_2  sum = CGAL::minkowski_sum_2(P, Q);
  std::cout << "P (+) Q = ";
  print_polygon_with_holes(sum);
  return 0;
}
