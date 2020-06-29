
#ifndef OCTREE_IO_H
#define OCTREE_IO_H

#include <CGAL/Octree.h>
#include <CGAL/Octree/Tree_walker_criterion.h>

#include <iostream>
#include <ostream>

using std::ostream;

template<class Kernel, class PointRange>
ostream &operator<<(ostream &os, const CGAL::Octree_node<Kernel, PointRange> &node) {

  // Show the depth of the node
  for (int i = 0; i < node.depth(); ++i)
    os << ". ";

  // Wrap information in brackets
  os << "{ ";

  // Index identifies which child this is
  os << "(" << node.index() << ") ";

  // Location
  os << "(" << node.location()[0] << ",";
  os << node.location()[1] << ",";
  os << node.location()[2] << ") ";

  //  os << "location: " << node.location() << ", ";
  if (!node.is_empty())
    os << "[" << node.num_points() << " points] ";

  os << (node.is_leaf() ? "[leaf] " : "");

  // Wrap information in brackets
  os << "}" << std::endl;

  return os;
}

template<class PointRange,
        class PointMap>
ostream &operator<<(ostream &os, const CGAL::Octree<PointRange, PointMap> &octree) {

  auto tree_walker = CGAL::Preorder();
  auto first = tree_walker.first(&octree.root());
  octree.print(os, first, tree_walker);
  return os;
}

#endif //OCTREE_IO_H
