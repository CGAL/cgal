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
#include <CGAL/boost/iterator/counting_iterator.hpp>

#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

template<typename InputIterator, typename PointMap>
struct Point_map_to_indexed_point_map {
  typedef std::size_t key_type;
  typedef typename boost::property_traits<PointMap>::value_type value_type;
  typedef typename boost::property_traits<PointMap>::reference reference;
  typedef typename boost::readable_property_map_tag category;

  InputIterator begin;
  PointMap point_map;

  Point_map_to_indexed_point_map(InputIterator begin = InputIterator(),
                                 PointMap point_map = PointMap())
          : begin(begin), point_map(point_map) {}

  friend reference get(const Point_map_to_indexed_point_map &map, std::size_t index) {
    return get(map.point_map, *(map.begin + index));
  }
};

template<class Traits>
class Direct_octree {

  typedef typename Traits::Input_range::iterator Input_iterator;
  typedef typename Traits::Point_map Point_map;
  typedef typename Traits::FT FT;
  typedef std::vector<std::size_t> Input_range;
  typedef Point_map_to_indexed_point_map<Input_iterator, Point_map> Indexed_point_map;

  typedef CGAL::Octree::Octree<Input_range, Indexed_point_map> Octree;

  Traits m_traits;
  Input_range m_input_range;
  Indexed_point_map m_index_map;
  CGAL::Octree::Octree<Input_range, Indexed_point_map> m_octree;

  std::size_t m_offset;

public:

  typedef typename Octree::Node Node;

  Direct_octree(const Traits &traits,
                Input_iterator begin,
                Input_iterator end,
                Point_map point_map,
                std::size_t offset = 0) :
          m_traits(traits),
          m_input_range(boost::counting_iterator<std::size_t>(0),
                        boost::counting_iterator<std::size_t>(end - begin)),
          m_index_map(begin, point_map),
          m_octree(m_input_range, m_index_map),
          m_offset(offset) {}

  std::size_t size() const {
    return m_octree.root().size();
  }

  std::size_t maxLevel() const {
    return m_octree.max_depth_reached();
  }

  std::size_t offset() const { return m_offset; }

  void refine(double cluster_epsilon_for_max_level_recomputation = -1., std::size_t bucketSize = 2,
              std::size_t maxLevel = 10) {

    if (cluster_epsilon_for_max_level_recomputation > 0.) {

      auto m_bBox = m_octree.bbox(m_octree.root());

      FT bbox_diagonal = (FT) CGAL::sqrt(
              (m_bBox.xmax() - m_bBox.xmin()) * (m_bBox.xmax() - m_bBox.xmin())
              + (m_bBox.ymax() - m_bBox.ymin()) * (m_bBox.ymax() - m_bBox.ymin())
              + (m_bBox.zmax() - m_bBox.zmin()) * (m_bBox.zmax() - m_bBox.zmin()));

      maxLevel = std::size_t(std::log(bbox_diagonal
                                      / cluster_epsilon_for_max_level_recomputation)
                             / std::log(2.0));

    }

    m_octree.refine(maxLevel, bucketSize);
  }

  typename Traits::FT width() const {
    return m_octree.bbox(m_octree.root()).xmax() - m_octree.bbox(m_octree.root()).xmin();
  }

  const Node &locate(const typename Traits::Point_3 &p) const {
    return m_octree.locate(p);
  }

  const Node &root() const { return m_octree.root(); }

  typename Traits::Point_3 barycenter(const Node &node) const {
    return m_octree.barycenter(node);
  }
};

template<class Traits>
class Indexed_octree {

  typedef typename Traits::Input_range::iterator Input_iterator;
  typedef typename Traits::Point_map Point_map;
  typedef std::vector<std::size_t> Input_range;
  typedef Point_map_to_indexed_point_map<Input_iterator, Point_map> Indexed_point_map;

  typedef Octree::Octree<Input_range, Indexed_point_map> Octree;

  Traits m_traits;
  Input_range m_input_range;
  Indexed_point_map m_index_map;
  Octree m_octree;

public:

  typedef typename Octree::Node Node;

  Indexed_octree(const Traits &traits,
                 Input_iterator begin,
                 Input_iterator end,
                 Point_map point_map,
                 std::size_t offset = 0) :
          m_traits(traits),
          m_input_range(boost::counting_iterator<std::size_t>(0),
                        boost::counting_iterator<std::size_t>(end - begin)),
          m_index_map(begin, point_map),
          m_octree(m_input_range, m_index_map, 1.0) {}

  std::size_t size() const {
    return m_octree.root().size();
  }

  Bbox_3 bbox() const {
    return m_octree.bbox(m_octree.root());
  }

  std::size_t maxLevel() const {
    return m_octree.max_depth_reached();
  }

  std::size_t offset() const { return 0; }

  void refine(double cluster_epsilon_for_max_level_recomputation = -1., std::size_t bucketSize = 2,
              std::size_t maxLevel = 10) {

    // TODO: I need to find out what cluster_epsilon is used for
    m_octree.refine(maxLevel, bucketSize);
  }

  std::size_t index(std::size_t i) { return m_index_map[i]; }

  typename Traits::FT width() const {
    return m_octree.bbox(m_octree.root()).xmax() - m_octree.bbox(m_octree.root()).xmin();
  }

  const Node &locate(const typename Traits::Point_3 &p) const {
    return m_octree.locate(p);
  }

  const Node &root() const { return m_octree.root(); }

  typename Traits::Point_3 barycenter(const Node &node) const {
    return m_octree.barycenter(node);
  }
};

}
}
}

#endif
