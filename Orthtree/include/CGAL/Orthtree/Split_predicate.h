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

#ifndef CGAL_ORTHTREE_SPLIT_PREDICATE_H
#define CGAL_ORTHTREE_SPLIT_PREDICATE_H

#include <CGAL/license/Orthtree.h>

#include <iostream>

namespace CGAL {

namespace Orthtrees {

namespace Split_predicate {

/*!
  \ingroup PkgOrthtreeSplitPredicates
  \brief predicate to split nodes of an orthtree when they contain more than a certain number of items
 */
class Bucket_size {

  std::size_t m_bucket_size;

public:
  /*!
    \brief creates a bucket size predicate.
   */
  Bucket_size(std::size_t bucket_size) :
          m_bucket_size(bucket_size) {}

  /*!
    \brief returns `true` if `n` should be splitted, `false` otherwise.
   */
  template<class Node>
  bool operator()(const Node &n) const {
    return (n.size() > m_bucket_size);
  }
};

/*!
  \ingroup PkgOrthtreeSplitPredicates
  \brief predicate to split nodes of an orthtree when they are less than a certain depth.
 */
class Max_depth {

  std::size_t m_max_depth;

public:

  /*!
    \brief creates a max depth predicate.
   */
  Max_depth(std::size_t max_depth) : m_max_depth(max_depth) {}

  /*!
    \brief returns `true` if `n` should be splitted, `false` otherwise.
   */
  template<class Node>
  bool operator()(const Node &n) const {
    return n.depth() < m_max_depth;
  }
};

/*!
  \ingroup PkgOrthtreeSplitPredicates
  \brief predicate to split nodes when they are less than a depth and they contain more than a number of items.
 */
class Max_depth_or_bucket_size {

  std::size_t m_max_depth, m_bucket_size;

public:

  /*!
    \brief creates a predicate using max depth or bucket size.
   */
  Max_depth_or_bucket_size(std::size_t max_depth, std::size_t bucket_size) :
          m_max_depth(max_depth), m_bucket_size(bucket_size) {}

  /*!
    \brief returns `true` if `n` should be splitted, `false` otherwise.
   */
  template<class Node>
  bool operator()(const Node &n) const {
    size_t num_points = n.size();
    size_t depth = n.depth();
    return (num_points > m_bucket_size && depth < m_max_depth);
  }
};

} // Split_predicate
} // Orthtrees
} // CGAL

#endif //CGAL_ORTHTREE_SPLIT_PREDICATE_H
