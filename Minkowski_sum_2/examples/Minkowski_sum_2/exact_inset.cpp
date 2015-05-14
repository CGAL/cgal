//! \file examples/Minkowski_sum_2/exact_inset.cpp
// Computing the exact inner offset of a polygon.

#include <iostream>

#include <CGAL/basic.h>

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
typedef Gps_traits::Polygon_2                   Offset_polygon;

int main(int argc, char* argv[])
{
  // Open the input file and read the input polygon.
  const char* filename = (argc > 1) ? argv[1] : "tight.dat";
  std::ifstream in_file(filename);
  if (! in_file.is_open()) {
    std::cerr << "Failed to open the input file." << std::endl;
    return -1;
  }
  Polygon_2 P;
  in_file >> P;
  in_file.close();
  std::cout << "Read an input polygon with " << P.size() << " vertices."
            << std::endl;

  // Compute the inner offset of the polygon.
  Traits traits;
  std::list<Offset_polygon> inset_polygons;
  boost::timer timer;
  inset_polygon_2(P, 1, traits, std::back_inserter(inset_polygons));
  double secs = timer.elapsed();

  std::list<Offset_polygon>::iterator it;
  std::cout << "The inset comprises "
            << inset_polygons.size() << " polygon(s)." << std::endl;
  for (it = inset_polygons.begin(); it != inset_polygons.end(); ++it)
      std::cout << "    Polygon with " << it->size() << " vertices."
                << std::endl;
  std::cout << "Inset computation took " << secs << " seconds." << std::endl;
  return 0;
}

#endif
