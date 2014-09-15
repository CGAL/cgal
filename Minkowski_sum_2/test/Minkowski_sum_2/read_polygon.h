#ifndef CGAL_READ_POLYGON_H
#define CGAL_READ_POLYGON_H

#include <CGAL/Polygon_2.h>
#include <iostream>
#include <fstream>

template <class Kernel>
void read_polygon (const char *filename, CGAL::Polygon_2<Kernel>& pgn)
{
  std::ifstream file(filename);

  if (!file)
  {
    std::cerr << "Failed to open " << filename << std::endl;
    exit(1);
  }

  file >> pgn;
}

#endif
