// Copyright (c) 2007-2020  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jackson Campolattaro, Cédric Portaneri, Tong Zhao

#ifndef CGAL_OCTREE_IO_H
#define CGAL_OCTREE_IO_H

#include <CGAL/license/Octree.h>

#include <CGAL/Octree.h>
#include <CGAL/Octree/Node.h>
#include <CGAL/Octree/Walker_criterion.h>

#include <iostream>
#include <ostream>
#include "Node.h"

using std::ostream;


template<typename Value>
ostream &operator<<(ostream &os, const CGAL::Octree::Node::Node<Value> &node) {

  // Show the depth of the node
//  for (int i = 0; i < node.depth(); ++i)
//    os << ". ";

  // Wrap information in brackets
  os << "{ ";

  // Index identifies which child this is
  os << "("
     << node.index()[0]
     << node.index()[1]
     << node.index()[2]
     << ") ";

  // Location
  os << "("
     << node.location()[0] << ","
     << node.location()[1] << ","
     << node.location()[2]
     << ") ";

  // Depth
  os << "("
     << +node.depth() // The + forces printing as an int instead of a char
     << ") ";

//  // If a node has points, indicate how many
//  if (!node.is_empty())
//    os << "[" << node.num_points() << " points] ";

  // If a node is a leaf, mark it
  os << (node.is_leaf() ? "[leaf] " : "");

  // If a node is root, mark it
  os << (node.is_root() ? "[root] " : "");

  // Wrap information in brackets
  os << "}";

  return os;
}

template<class PointRange,
        class PointMap>
ostream &operator<<(ostream &os, const CGAL::Octree::Octree<PointRange, PointMap> &octree) {

  // Create a range of nodes
//  auto nodes = octree.template walk<CGAL::Octree::Walker::Preorder>();
  auto nodes = octree.walk(CGAL::Octree::Walker::Preorder());

  // Iterate over the range and print each node
//  for (auto &n : nodes) {
//
//    for (int i = 0; i < n.depth() - 1; ++i)
//      os << " │  ";
//
//    if (!n.is_root()) {
//
//      if (n.index() == 7)
//        os << " └─";
//      else
//        os << " ├─";
//    }
//
//    os << n << std::endl;
//
//    if (!n.is_leaf()) {
//
//      for (int i = 0; i < n.depth(); ++i)
//        os << " │  ";
//
//      os << " ┬ " << std::endl;
//    }
//  }

  // Iterate over the range
  for (auto &n : nodes) {

    // Show the depth
    for (int i = 0; i < n.depth(); ++i)
      os << ". ";

    // Print the node
    os << n << std::endl;
  }

  return os;
}

#endif //CGAL_OCTREE_IO_H
