// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_OCTREE_H
#define CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_OCTREE_H

#include <CGAL/license/Shape_detection.h>

#include <stack>
#include <limits>

#include <CGAL/Random.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Shape_base.h>

#include <CGAL/Octree.h>
#include <CGAL/Octree/Node.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {


template<class Traits>
class Direct_octree {

  typedef typename Traits::Input_range::iterator Input_iterator;
  typedef typename Traits::Point_map Point_map;
  typedef std::vector<std::size_t> Input_range;

  typedef Octree::Octree<Input_range, typename Traits::Point_map> Octree;

  Traits m_traits;
  std::size_t m_offset;

  Octree m_octree;

public:

  typedef typename Octree::Node Node;

  Direct_octree(const Traits &traits,
                const Input_iterator &begin,
                const Input_iterator &end,
                Point_map &point_map,
                std::size_t offset = 0) :
          m_octree({begin, end}, point_map),
          m_traits(traits),
          m_offset(offset) {

  }

  std::size_t size() const {
    return m_octree.root().size();
  }

  const Bbox_3 &boundingBox() const {
    return m_octree.bbox(this->root());
  }

  std::size_t maxLevel() const {
    return m_octree.max_depth_reached();
  }

  std::size_t offset() const { return m_offset; }

  void createTree(double cluster_epsilon_for_max_level_recomputation = -1., std::size_t bucketSize = 2,
                  std::size_t maxLevel = 10) {

    // TODO: I need to find out what cluster_epsilon is used for
    m_octree.refine(maxLevel, bucketSize);
  }

  typename Traits::FT width() const { return m_octree.m_side_per_depth[0]; }

  const Node &locate(const typename Traits::Point_3 &p) const {
    return m_octree.locate(p);
  }

  const Node &root() const { return m_octree.root(); }
};

template<class Traits>
class Indexed_octree {

  typedef typename Traits::Input_range::iterator Input_iterator;
  typedef typename Traits::Point_map Point_map;
  typedef std::vector<std::size_t> Input_range;

  typedef Octree::Octree<Input_range, typename Traits::Point_map> Octree;

  Traits m_traits;
  std::vector<std::size_t> m_index_map;

  Octree m_octree;

public:

  typedef typename Octree::Node Node;

  Indexed_octree(const Traits &traits,
                 const Input_iterator &begin,
                 const Input_iterator &end,
                 Point_map &point_map) :
          m_octree({}, point_map),
          m_traits(traits) {

  }

  std::size_t size() const {
    return m_octree.root().size();
  }

  const Bbox_3 &boundingBox() const {
    return m_octree.bbox(m_octree.root());
  }

  std::size_t maxLevel() const {
    return m_octree.max_depth_reached();
  }

  void createTree(double cluster_epsilon_for_max_level_recomputation = -1., std::size_t bucketSize = 2,
                  std::size_t maxLevel = 10) {

    // TODO: I need to find out what cluster_epsilon is used for
    m_octree.refine(maxLevel, bucketSize);
  }

  std::size_t index(std::size_t i) { return m_index_map[i]; }

  typename Traits::FT width() const { return m_octree.m_side_per_depth[0]; }

  const Node &locate(const typename Traits::Point_3 &p) const {
    return m_octree.locate(p);
  }

  const Node &root() const { return m_octree.root(); }
};

}
}
}

#endif
