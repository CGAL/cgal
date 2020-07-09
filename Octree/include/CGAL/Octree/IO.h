
#ifndef OCTREE_IO_H
#define OCTREE_IO_H

#include <CGAL/Octree.h>
#include <CGAL/Octree/Node.h>
#include <CGAL/Octree/Walker_criterion.h>

#include <iostream>
#include <ostream>
#include "Node.h"

using std::ostream;

template<class PointRange,
        class PointMap>
ostream &operator<<(ostream &os, const CGAL::Octree::Octree<PointRange, PointMap> &octree) {

  // Create a range of nodes
  auto nodes = octree.template walk<CGAL::Octree::Walker::Preorder_tree_walker>();

  // Iterate over the range and print each node
  for (auto &n : nodes)
    os << n;

  return os;
}

template<typename Value>
ostream &operator<<(ostream &os, const CGAL::Octree::Node::Node<Value> &node) {

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

//  // If a node has points, indicate how many
//  if (!node.is_empty())
//    os << "[" << node.num_points() << " points] ";

  // If a node is a leaf, mark it
  os << (node.is_leaf() ? "[leaf] " : "");

  // If a node is root, mark it
  os << (node.is_root() ? "[root] " : "");

  // Wrap information in brackets
  os << "}" << std::endl;

  return os;
}

#endif //OCTREE_IO_H
