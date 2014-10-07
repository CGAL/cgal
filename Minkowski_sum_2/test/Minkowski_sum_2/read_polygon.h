#ifndef CGAL_READ_POLYGON_H
#define CGAL_READ_POLYGON_H

#include <CGAL/Polygon_with_holes_2.h>
#include <iostream>
#include <fstream>

template <typename Kernel>
bool read_polygon(const char* filename, CGAL::Polygon_2<Kernel>& pgn)
{
  std::ifstream file(filename);
  if (!file) {
    std::cerr << "Failed to open " << filename << std::endl;
    return false;
  }

  file >> pgn;
  return true;
}

template <typename Kernel>
bool read_polygon(const char* filename, CGAL::Polygon_with_holes_2<Kernel>& pgn)
{
  std::ifstream file(filename);
  if (!file) {
    std::cerr << "Failed to open " << filename << std::endl;
    return false;
  }

  file >> pgn;
  return true;
}

#endif
