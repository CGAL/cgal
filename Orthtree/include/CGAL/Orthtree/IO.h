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

#ifndef CGAL_ORTHTREE_IO_H
#define CGAL_ORTHTREE_IO_H

#include <CGAL/license/Orthtree.h>

#include <iostream>
#include <ostream>

namespace CGAL
{
namespace internal
{

template<typename Node_index, typename Tree>
std::ostream& print_orthtree_node(std::ostream& os, const Node_index& node, const Tree &tree)
{

  // Wrap information in brackets
  os << "{ ";

  // Index identifies which child this is
  os << "(";
  for (std::size_t i = 0; i < tree.local_coordinates(node).size(); ++ i)
    os << tree.local_coordinates(node)[i];
  os << ") ";

  // Location
  os << "( ";
  for (const auto& b : tree.global_coordinates(node))
    os << b << " ";
  os << ") ";

  // Depth
  os << "("
     << +tree.depth(node) // The + forces printing as an int instead of a char
     << ") ";

  os << "("
     << tree.data(node).size()
     << ") ";

  // If a node is a leaf, mark it
  os << (tree.is_leaf(node) ? "[leaf] " : "");

  // If a node is root, mark it
  os << (tree.is_root(node) ? "[root] " : "");

  // Wrap information in brackets
  os << "}";

  return os;
}

} // internal
} // CGAL

#endif //CGAL_ORTHTREE_IO_H
