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
#include <CGAL/Shape_detection/Efficient_RANSAC/Efficient_RANSAC_traits.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>

#include <CGAL/Octree.h>

namespace CGAL {
namespace Shape_detection {

// Forward declaration needed for automatic traits detection without
// including the deprecated header itself…
template <typename Gt, typename IR, typename IPM, typename INM>
struct Shape_detection_traits;

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

template <typename Traits>
struct Traits_base { typedef Traits type; };
template <typename Gt, typename IR, typename IPM, typename INM>
struct Traits_base<CGAL::Shape_detection::Efficient_RANSAC_traits<Gt,IR,IPM,INM> >
{ typedef Gt type; };
template <typename Gt, typename IR, typename IPM, typename INM>
struct Traits_base<CGAL::Shape_detection::Shape_detection_traits<Gt,IR,IPM,INM> >
{ typedef Gt type; };

template<class Traits>
class RANSAC_octree {

  typedef typename Traits::Input_range::iterator Input_iterator;
  typedef typename Traits::Point_map Point_map;
  typedef typename Traits::FT FT;
  typedef std::vector<std::size_t> Input_range;
  typedef Point_map_to_indexed_point_map<Input_iterator, Point_map> Indexed_point_map;

  typedef CGAL::Octree<typename Traits_base<Traits>::type,
                       Input_range, Indexed_point_map> Octree;

  Traits m_traits;
  Input_range m_input_range;
  Indexed_point_map m_index_map;
  Octree m_octree;

  std::size_t m_offset;

public:

  typedef typename Octree::Node Node;

  RANSAC_octree(const Traits &traits,
                Input_iterator begin,
                Input_iterator end,
                Point_map point_map,
                std::size_t offset = 0) :
          m_traits(traits),
          m_input_range(boost::counting_iterator<std::size_t>(0),
                        boost::counting_iterator<std::size_t>(end - begin)),
          m_index_map(begin, point_map),
          m_octree(m_input_range, m_index_map, 1.0),
          m_offset(offset) {}

  std::size_t index (Node node, std::size_t i) const
  {
    return m_offset + *(node.points().begin() + i);
  }

  std::size_t size() const {
    return m_octree.root().size();
  }

  std::size_t maxLevel() const {
    return m_octree.depth() - 1;
  }

  std::size_t offset() const { return m_offset; }

  void refine(double cluster_epsilon_for_max_level_recomputation = -1.,
              std::size_t bucketSize = 20,
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

  Node locate(const typename Traits::Point_3 &p) const {
    return m_octree.locate(p);
  }

  Node root() const { return m_octree.root(); }

  typename Traits::Point_3 barycenter(const Node &node) const {
    return m_octree.barycenter(node);
  }

  Bbox_3 boundingBox() const {
      return m_octree.bbox(m_octree.root());
  }
};

}
}
}

#endif
