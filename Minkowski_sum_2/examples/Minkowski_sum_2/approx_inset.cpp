//! \file examples/Minkowski_sum_2/approx_inset.cpp
// Computing the approximated inset of a polygon.

#include <fstream>
#include <iostream>
#include <list>
#include <boost/timer.hpp>

#include <CGAL/basic.h>
#include <CGAL/approximated_offset_2.h>

#include "bops_circular.h"

typedef CGAL::Polygon_2<Kernel>                         Linear_polygon;

int main(int argc, char* argv[])
{
  // Open the input file and read a polygon.
  const char* filename = (argc > 1) ? argv[1] : "tight.dat";
  std::ifstream in_file(filename);

  if (! in_file.is_open()) {
    std::cerr << "Failed to open the input file." << std::endl;
    return -1;
  }

  // Read the input polygon.
  Linear_polygon P;
  in_file >> P;
  in_file.close();

  std::cout << "Read an input polygon with " << P.size() << " vertices."
            << std::endl;

  // Approximate the offset polygon.
  std::list<Polygon_2> inset_polygons;
  boost::timer timer;
  approximated_inset_2(P, 1, 0.00001, std::back_inserter(inset_polygons));
  double secs = timer.elapsed();

  std::list<Polygon_2>::iterator it;
  std::cout << "The inset comprises " << inset_polygons.size()
            << " polygon(s)." << std::endl;
  for (it = inset_polygons.begin(); it != inset_polygons.end(); ++it)
    std::cout << "    Polygon with " << it->size() << " vertices." << std::endl;
  std::cout << "Inset computation took " << secs << " seconds." << std::endl;
  return 0;
}
