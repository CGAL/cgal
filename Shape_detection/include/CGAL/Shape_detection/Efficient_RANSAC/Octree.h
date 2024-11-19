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
#include <CGAL/property_map.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>

#include <CGAL/Octree.h>

namespace CGAL {
namespace Shape_detection {


namespace internal {

template <typename Traits>
struct Traits_base { typedef Traits type; };
template <typename Gt, typename IR, typename IPM, typename INM>
struct Traits_base<CGAL::Shape_detection::Efficient_RANSAC_traits<Gt,IR,IPM,INM> >
{ typedef Gt type; };

template<class Traits>
class RANSAC_octree {

  typedef typename Traits::Input_range::iterator Input_iterator;
  typedef typename Traits::Point_map Point_map;
  typedef typename Traits::FT FT;
  typedef std::vector<std::size_t> Input_range;
  typedef Random_index_access_property_map<Input_iterator, Point_map> Indexed_point_map;

  typedef Orthtree_traits_point<typename Traits_base<Traits>::type, Input_range, Indexed_point_map, false, 3> OTraits;

  typedef CGAL::Orthtree<OTraits> Octree;

  Traits m_traits;
  Input_range m_input_range;
  Indexed_point_map m_index_map;
  Octree m_octree;
  Bbox_3 m_bBox;
  FT m_width;

  std::size_t m_offset;

public:

  typedef typename Octree::Node_index Node;
  typedef typename OTraits::Node_data Node_data;

  RANSAC_octree(const Traits &traits,
                Input_iterator begin,
                Input_iterator end,
                Point_map point_map,
                std::size_t offset = 0) :
          m_traits(traits),
          m_input_range(boost::counting_iterator<std::size_t>(0),
                        boost::counting_iterator<std::size_t>(end - begin)),
          m_index_map(begin, point_map),
          m_octree(OTraits(m_input_range, m_index_map)),
          m_bBox (bbox_3(make_transform_iterator_from_property_map(begin, point_map),
                         make_transform_iterator_from_property_map(end, point_map))),
          m_offset(offset) {}

  std::size_t index (Node node, std::size_t i) const
  {
    return m_offset + *(m_octree.data(node).begin() + i);
  }

  std::size_t depth(const Node& node) const {
    return m_octree.depth(node);
  }

  bool is_leaf(const Node& node) const {
    return m_octree.is_leaf(node);
  }

  std::size_t size() const {
    return m_input_range.size();
  }

  std::size_t maxLevel() const {
    return m_octree.depth();
  }

  std::size_t offset() const { return m_offset; }

  void refine(double cluster_epsilon_for_max_level_recomputation = -1.,
              std::size_t bucketSize = 20,
              std::size_t maxLevel = 10) {

    if (cluster_epsilon_for_max_level_recomputation > 0.) {

      FT bbox_diagonal = (FT) CGAL::sqrt(
              (m_bBox.xmax() - m_bBox.xmin()) * (m_bBox.xmax() - m_bBox.xmin())
              + (m_bBox.ymax() - m_bBox.ymin()) * (m_bBox.ymax() - m_bBox.ymin())
              + (m_bBox.zmax() - m_bBox.zmin()) * (m_bBox.zmax() - m_bBox.zmin()));

      maxLevel = std::size_t(std::log2(bbox_diagonal
                                       / cluster_epsilon_for_max_level_recomputation));
      if (maxLevel == 0)
        maxLevel = 1;
    }

    m_octree.refine(maxLevel, bucketSize);

    m_width = FT(0.5) * FT(m_octree.bbox(m_octree.root()).xmax() - m_octree.bbox(m_octree.root()).xmin());
  }

  const typename Traits::FT& width() const {
    return m_width;
  }

  Node child(const Node& node, std::size_t i) const {
    return m_octree.child(node, i);
  }

  Node parent(const Node& node) const {
    return m_octree.parent(node);
  }

  Node locate(const typename Traits::Point_3 &p) const {
    return m_octree.locate(p);
  }

  Node root() const { return m_octree.root(); }

  Node_data points(const Node& n) const { return m_octree.data(n); }

  typename Traits::Point_3 barycenter(const Node &node) const {
    return m_octree.barycenter(node);
  }

  typename Traits::GeomTraits::Iso_cuboid_3 boundingBox() const {
      return m_octree.bbox(m_octree.root());
  }
};

}
}
}

#endif
