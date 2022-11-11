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

#ifndef CGAL_ORTHTREE_SPLIT_PREDICATES_H
#define CGAL_ORTHTREE_SPLIT_PREDICATES_H

#include <CGAL/license/Orthtree.h>

#include <iostream>

namespace CGAL {

namespace Orthtrees {

/*!
  \ingroup PkgOrthtreeSplitPredicates
  \brief A class used to choose when a node should be split depending on the number of inliers.

  This is a bucket size predicate that considers a node should be
  split if it contains more than a certain number of items.
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
    \brief returns `true` if `n` should be split, `false` otherwise.
   */
  template<typename Node>
  bool operator()(const Node &n) const {
    return (n.size() > m_bucket_size);
  }
};

/*!
  \ingroup PkgOrthtreeSplitPredicates
  \brief A class used to choose when a node should be split depending on the depth.

  This predicate makes a node be split if its depth is lower than a certain limit.
 */
class Maximum_depth {

  std::size_t m_max_depth;

public:

  /*!
    \brief creates a maximum depth predicate.
   */
  Maximum_depth(std::size_t max_depth) : m_max_depth(max_depth) {}

  /*!
    \brief returns `true` if `n` should be split, `false` otherwise.
   */
  template<typename Node>
  bool operator()(const Node &n) const {
    return n.depth() < m_max_depth;
  }
};

/*!
  \ingroup PkgOrthtreeSplitPredicates
  \brief A class used to choose when a node should be split depending on the depth and the number of inliers.

  This predicate makes a note split if it contains more than a
  certain number of items and if its depth is lower than a certain
  limit.

  The refinement is stopped as soon as one of the conditions is
  violated: if a node has more inliers than `bucket_size` but is
  already at `max_depth`, it is not split. Similarly, a node that is
  at a depth smaller than `max_depth` but already has fewer inliers
  than `bucket_size`, it is not split.

 */
class Maximum_depth_and_maximum_number_of_inliers {

  std::size_t m_max_depth, m_bucket_size;

public:

  /*!  \brief creates a predicate using maximum depth or bucket size.
   */
  Maximum_depth_and_maximum_number_of_inliers(std::size_t max_depth, std::size_t bucket_size) :
          m_max_depth(max_depth), m_bucket_size(bucket_size) {}

  /*!
    \brief returns `true` if `n` should be split, `false` otherwise.
   */
  template<typename Node>
  bool operator()(const Node &n) const {
    std::size_t num_points = n.size();
    std::size_t depth = n.depth();
    return (num_points > m_bucket_size && depth < m_max_depth);
  }
};

} // Orthtrees
} // CGAL

#endif //CGAL_ORTHTREE_SPLIT_PREDICATES_H
