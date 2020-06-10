
#ifndef OCTREE_IO_H
#define OCTREE_IO_H

#include <CGAL/Octree.h>

#include <iostream>

using std::ostream;

template<class Kernel, class PointRange>
ostream& operator<<(ostream& os, const CGAL::Octree_node<Kernel, PointRange>& node) {

  os << "[ depth: " << node.depth() << " ]";

  return os;
}

#endif //OCTREE_IO_H
