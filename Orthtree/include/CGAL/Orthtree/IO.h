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

template<typename Node>
std::ostream& print_orthtree_node(std::ostream& os, const Node& node)
{
  // Show the depth of the node
//  for (int i = 0; i < node.depth(); ++i)
//    os << ". ";

  // Wrap information in brackets
  os << "{ ";

  // Index identifies which child this is
  os << "(";
  for (std::size_t i = 0; i < node.local_coordinates().size(); ++ i)
    os << node.local_coordinates()[i];
  os << ") ";

  // Location
  os << "( ";
  for (const auto& b : node.global_coordinates())
    os << b << " ";
  os << ") ";

  // Depth
  os << "("
     << +node.depth() // The + forces printing as an int instead of a char
     << ") ";

  os << "("
     << node.size()
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

} // internal
} // CGAL

#endif //CGAL_ORTHTREE_IO_H
