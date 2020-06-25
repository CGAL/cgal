
#ifndef OCTREE_IO_H
#define OCTREE_IO_H

#include <CGAL/Octree.h>

#include <iostream>
#include <ostream>

using std::ostream;

template<class Kernel, class PointRange>
void writeToStream(ostream &os, const CGAL::Octree_node<Kernel, PointRange> &node) {

  // Set indentation
  for (int i = 0; i < node.depth(); ++i) {
    os << ". ";
  }

  // Print out this node
  os << "(" << node.location()[0] << "," << node.location()[1] << "," << node.location()[2] << ")";

  if (node.num_points() > 0)
    os << " [" << node.num_points() << " points]";

  os << std::endl;

  // Print out this node's children
  if (!node.is_leaf()) {
    for (int i = 0; i < 8; ++i) {
      writeToStream(os, node[i]);
    }
  }
}

template<class Kernel, class PointRange>
ostream &operator<<(ostream &os, const CGAL::Octree_node<Kernel, PointRange> &node) {

  writeToStream(os, node);
  return os;
}

#endif //OCTREE_IO_H
