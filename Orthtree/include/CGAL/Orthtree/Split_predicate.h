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
  \brief bucket size predicate that splits a node if it contains more than a certain number of items.
 */
class Maximum_number_of_inliers {

  std::size_t m_bucket_size;

public:
  /*!
    \brief creates a predicate based on the number of inliers (bucket size).
   */
  Maximum_number_of_inliers(std::size_t bucket_size) :
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
  \brief predicate that splits a node if its depth is lower than a certain limit.
 */
class Maximum_depth {

  std::size_t m_max_depth;

public:

  /*!
    \brief creates a maximum depth predicate.
   */
  Maximum_depth(std::size_t max_depth) : m_max_depth(max_depth) {}

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

  \brief predicate that splits a node if it contains more than a
  certain number of items and if its depth is lower than a certain
  limit.
 */
class Maximum_depth_and_maximum_number_of_inliers {

  std::size_t m_max_depth, m_bucket_size;

public:

  /*!
    \brief creates a predicate using minimum depth or bucket size.
   */
  Maximum_depth_and_maximum_number_of_inliers(std::size_t max_depth, std::size_t bucket_size) :
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
