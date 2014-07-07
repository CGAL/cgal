#ifndef _OFF_H_
#define _OFF_H_

#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>

// Write into an OFF file to visualize with GeomView
template <typename K, typename Polyhedron>
void convertToOFF (std::string const& filename, Polyhedron& P) {
  std::ofstream file(filename.c_str());
  file << P;
}

#endif

