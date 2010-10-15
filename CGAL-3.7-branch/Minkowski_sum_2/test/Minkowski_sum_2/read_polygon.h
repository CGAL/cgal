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
  // Open the input file.
  std::ifstream          ifile (filename);

  if (! ifile.is_open())
  {
    std::cerr << "Failed to open <" << filename << ">." << std::endl;
    return (false);
  }

  // Read the polygon.
  int                                     n_vertices = 0;
  typename Kernel::FT                     x, y;
  std::list<typename Kernel::Point_2>     vertices;
  int                                     k;

  // Read the number of polygon vertices.
  ifile >> n_vertices;

  // Read the vertices.
  for (k = 0; k < n_vertices; k++)
  {
    ifile >> x >> y;

    vertices.push_back (typename Kernel::Point_2 (x, y));
  }
  ifile.close();

  pgn = CGAL::Polygon_2<Kernel> (vertices.begin(), vertices.end());

  // Make sure the polygon is simple.
  if (! pgn.is_simple())
  {
    std::cerr << "Error - the polygon is not simple." << std::endl;
    return (false);
  }

  return (true);
}

#endif
