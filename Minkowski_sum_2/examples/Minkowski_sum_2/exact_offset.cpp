//! \file examples/Minkowski_sum_2/exact_offset.cpp
// Computing the exact offset of a polygon.

#include <iostream>

#include <CGAL/config.h>

#ifndef CGAL_USE_CORE
int main()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return 0;
}
#else

#include <fstream>
#include <boost/timer.hpp>

#include <CGAL/Gps_traits_2.h>
#include <CGAL/offset_polygon_2.h>

#include "arr_conics.h"

typedef CGAL::Polygon_2<Rat_kernel>             Polygon_2;
typedef CGAL::Gps_traits_2<Traits>              Gps_traits;
typedef Gps_traits::Polygon_with_holes_2        Offset_polygon_with_holes_2;

int main(int argc, char* argv[])
{
  // Open the input file and read the input polygon.
  const char* filename = (argc > 1) ? argv[1] : "spiked.dat";
  std::ifstream in_file(filename);
  if (! in_file.is_open()) {
    std::cerr << "Failed to open the input file." << std::endl;
    return -1;
  }
  Polygon_2  P;
  in_file >> P;
  in_file.close();
  std::cout << "Read an input polygon with " << P.size() << " vertices."
            << std::endl;

  // Compute the offset polygon.
  Traits traits;
  boost::timer timer;
  Offset_polygon_with_holes_2 offset = CGAL::offset_polygon_2(P, 5, traits);
  double secs = timer.elapsed();

  std::cout << "The offset polygon has " << offset.outer_boundary().size()
            << " vertices, " << offset.number_of_holes() << " holes."
            << std::endl;
  std::cout << "Offset computation took " << secs << " seconds." << std::endl;
  return 0;
}

#endif
