#ifndef CGAL_READ_POLYGON_TEST_H
#define CGAL_READ_POLYGON_TEST_H

#include <CGAL/Polygon_2.h>
#include <iostream>
#include <fstream>

/*!
 * Read a polygon from an input file.
 * \param filename The name of the input file.
 * \param pgn Output: The polygon.
 * \return Whether the polygon was successfuly read.
 */
template <class Kernel>
bool read_polygon (const char *filename, CGAL::Polygon_2<Kernel>& pgn)
{
  std::ifstream ifile(filename);
  ifile >> pgn;

  // Make sure the polygon is simple.
  if (! pgn.is_simple())
  {
    std::cerr << "Error - the polygon is not simple." << std::endl;
    return false;
  }

  return true;
}

template <class Kernel>
bool read_polygon_with_holes (const char *filename, CGAL::Polygon_with_holes_2<Kernel>& pgn)
{
  std::ifstream ifile(filename);
  ifile >> pgn;
  // TODO: what can go wrong?
  return true;
}

template <class Kernel>
bool write_polygon_with_holes (const char *filename, CGAL::Polygon_with_holes_2<Kernel>& pgn)
{
  std::cout << filename << std::endl;
  std::ofstream ofile(filename);
  if (ofile.is_open()) {
      ofile << pgn;
  } else {
      exit(1);
  }
  return true;
}

#endif
