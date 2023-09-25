// Copyright (c) 2023  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jackson Campolattaro

#ifndef ORTHTREE_TESTS_ORTHTREE_TRAITS_POINT_H
#define ORTHTREE_TESTS_ORTHTREE_TRAITS_POINT_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Dimension.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Point_set_2.h>
#include <CGAL/Orthtree/Cartesian_ranges.h>

#include <CGAL/Orthtree_traits_base_for_dimension.h>

namespace CGAL {

template <typename Tree, typename PointMap>
void reassign_points(
  Tree& tree, PointMap& point_map,
  typename Tree::Node_index n, const typename Tree::Point& center, typename Tree::Node_data points,
  std::bitset<Tree::Dimension::value> coord = {}, std::size_t dimension = 0
) {

  // Root case: reached the last dimension
  if (dimension == Tree::Dimension::value) {
    tree.data(tree.child(n, coord.to_ulong())) = points;
    return;
  }

  // Split the point collection around the center point on this dimension
  auto split_point = std::partition(
    points.begin(), points.end(),
    [&](const auto& p) -> bool {
      // This should be done with cartesian iterator,
      // but it seems complicated to do efficiently
      return (get(point_map, p)[int(dimension)] < center[int(dimension)]);
    }
  );

  // Further subdivide the first side of the split
  std::bitset<Tree::Dimension::value> coord_left = coord;
  coord_left[dimension] = false;
  reassign_points(tree, point_map, n, center, {points.begin(), split_point}, coord_left, dimension + 1);

  // Further subdivide the second side of the split
  std::bitset<Tree::Dimension::value> coord_right = coord;
  coord_right[dimension] = true;
  reassign_points(tree, point_map, n, center, {split_point, points.end()}, coord_right, dimension + 1);
}

/*!
  \ingroup PkgOrthtreeTraits

  The class `Orthtree_traits_point` can be used as a template parameter of
  the `Orthtree` class.

  \tparam GeomTraits model of `Kernel`.
  \tparam PointSet must be a model of range whose value type is the key type of `PointMap`
  \tparam PointMap must be a model of `ReadablePropertyMap` whose value type is `GeomTraits::Traits::Point_d`

  \warning The input point set is not copied. It is used directly
  and is rearranged by the `Orthtree`. Altering the point range
  after creating the orthtree might leave it in an invalid state.

  \cgalModels{OrthtreeTraits}
  \sa `CGAL::Octree`
  \sa `CGAL::Orthtree_traits_2`
  \sa `CGAL::Orthtree_traits_3`
  \sa `CGAL::Orthtree_traits_d`
*/
template <
  typename GeomTraits,
  typename PointSet,
  typename PointMap = Identity_property_map<typename std::iterator_traits<typename PointSet::iterator>::value_type>,
  typename DimensionTag = Ambient_dimension<
    typename std::iterator_traits<typename PointSet::iterator>::value_type,
    GeomTraits
  >
>
struct Orthtree_traits_point : public Orthtree_traits_base_for_dimension<GeomTraits, DimensionTag> {
public:

  /// \name Types
  /// @{

  using Self = Orthtree_traits_point<GeomTraits, PointSet, PointMap, DimensionTag>;
  using Tree = Orthtree<Self>;

  using Node_data = boost::iterator_range<typename PointSet::iterator>;
  using Node_data_element = typename std::iterator_traits<typename PointSet::iterator>::value_type;

  /// @}

  Orthtree_traits_point(
    PointSet& point_set,
    PointMap point_map = PointMap()
  ) : m_point_set(point_set), m_point_map(point_map) {}

  /// \name Operations
  /// @{

  auto construct_root_node_bbox_object() const {
    return [&]() -> typename Self::Bbox_d {

      std::array<typename Self::FT, Self::Dimension::value> bbox_min, bbox_max;
      Orthtrees::internal::Cartesian_ranges<Self> cartesian_range;

      // init bbox with first values found
      {
        const typename Self::Point_d& point = get(m_point_map, *(m_point_set.begin()));
        std::size_t i = 0;
        for (const typename Self::FT& x: cartesian_range(point)) {
          bbox_min[i] = x;
          bbox_max[i] = x;
          ++i;
        }
      }
      // Expand bbox to contain all points
      for (const auto& p: m_point_set) {
        const typename Self::Point_d& point = get(m_point_map, p);
        std::size_t i = 0;
        for (const typename Self::FT& x: cartesian_range(point)) {
          bbox_min[i] = (std::min)(x, bbox_min[i]);
          bbox_max[i] = (std::max)(x, bbox_max[i]);
          ++i;
        }
      }

      return {std::apply(Self::construct_point_d_object(), bbox_min),
              std::apply(Self::construct_point_d_object(), bbox_max)};
    };
  }

  auto construct_root_node_contents_object() const {
    return [&]() -> typename Self::Node_data {
      return {m_point_set.begin(), m_point_set.end()};
    };
  }

  auto distribute_node_contents_object() const {
    return [&](typename Tree::Node_index n, Tree& tree, const typename Self::Point_d& center) {
      CGAL_precondition(!tree.is_leaf(n));
      reassign_points(tree, m_point_map, n, center, tree.data(n));
    };
  }

  auto get_geometric_object_for_element_object() const {
    return [&](const Node_data_element& index) -> typename Self::Point_d {
      return get(m_point_map, index);
    };
  }

  /// @}

private:

  PointSet& m_point_set;
  PointMap m_point_map;

};

}


#endif //ORTHTREE_TESTS_ORTHTREE_TRAITS_POINT_H
