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

#ifndef CGAL_OCTREE_SPLIT_CRITERION_H
#define CGAL_OCTREE_SPLIT_CRITERION_H

#include <CGAL/license/Octree.h>

#include <iostream>

namespace CGAL {

namespace Octree {

namespace Split_criterion {

/*!
 * \brief Criterion to split nodes of an octree when they contain more than a certain number of items
 */
struct Bucket_size {

  std::size_t m_bucket_size;

  Bucket_size(std::size_t bucket_size) :
          m_bucket_size(bucket_size) {}

  template<class Node>
  bool operator()(const Node &n) const {
    return (n.number_of_points() > m_bucket_size);
  }
};

/*!
 * \brief Criterion to split nodes of an octree when they are less than a certain depth
 */
struct Max_depth {

  std::size_t m_max_depth;

  Max_depth(std::size_t max_depth) : m_max_depth(max_depth) {}

  template<class Node>
  bool operator()(const Node &n) const {
    return n.depth() < m_max_depth;
  }
};

/*!
 * \brief Criterion to split nodes when they are less than a depth and they contain more than a number of items
 */
struct Max_depth_or_bucket_size {

  std::size_t m_max_depth, m_bucket_size;

  Max_depth_or_bucket_size(std::size_t max_depth, std::size_t bucket_size) :
          m_max_depth(max_depth), m_bucket_size(bucket_size) {}

  template<class Node>
  bool operator()(const Node &n) const {
    size_t num_points = n.number_of_points();
    size_t depth = n.depth();
    return (num_points > m_bucket_size && depth < m_max_depth);
  }
};

}
}
}

/*

// Just for an example using outside info (like a normal map)
// The normals remain unknown to the octree but are used for construction
struct Stop_at_normal_deviation {
 Normal_map normal_map;
 FT max_dev;

 Stop_at_normal_deviation(Normal_map normal_map,
                          FT max_dev)
         : normal_map(normal_map), max_dev(max_dev) {}

 bool operator()(const Node &n) const {
   FT dev = 0;

   for (Iterator it = n.begin(); it != n.end(); ++it)
     dev += compute_deviation(get(normal_map, *it)); // whatever compute_deviation is :)

   // if your node defines begin() and end(), you can also use a C++11 loop:
   // for (const Range_type& r : n)
   //   dev += compute_deviation(get (normal_map, r));

   return dev < max_dev;
 }
};

*/

#endif //CGAL_OCTREE_SPLIT_CRITERION_H
