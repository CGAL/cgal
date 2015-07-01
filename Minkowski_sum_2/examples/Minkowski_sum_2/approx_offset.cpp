//! \file examples/Minkowski_sum_2/approx_offset.cpp
// Computing the approximated offset of a polygon.

#include <fstream>
#include <boost/timer.hpp>

#include <CGAL/basic.h>
#include <CGAL/approximated_offset_2.h>

#include "bops_circular.h"

typedef CGAL::Polygon_2<Kernel>                         Linear_polygon;

int main(int argc, char* argv[])
{
  // Open the input file and read a polygon.
  const char* filename = (argc > 1) ? argv[1] : "spiked.dat";
  std::ifstream in_file(filename);
  if (! in_file.is_open()) {
    std::cerr << "Failed to open the input file." << std::endl;
    return -1;
  }
  Linear_polygon  P;
  in_file >> P;
  in_file.close();
  std::cout << "Read an input polygon with " << P.size() << " vertices."
            << std::endl;

  // Approximate the offset polygon with radius 5 and error bound 0.00001.
  boost::timer timer;
  Polygon_with_holes_2 offset = CGAL::approximated_offset_2(P, 5, 0.00001);
  double secs = timer.elapsed();

  std::cout << "The offset polygon has " << offset.outer_boundary().size()
            << " vertices, " << offset.number_of_holes() << " holes."
            << std::endl;
  std::cout << "Offset computation took " << secs << " seconds." << std::endl;
  return 0;
}
