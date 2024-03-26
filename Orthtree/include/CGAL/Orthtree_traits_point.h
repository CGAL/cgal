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

#include <CGAL/Orthtree_traits_base.h>

namespace CGAL {

template <typename Tree, typename PointMap>
void reassign_points(
  Tree& tree, PointMap& point_map,
  typename Tree::Node_index n, const typename Tree::Point& center, typename Tree::Node_data points,
  std::bitset<Tree::dimension> coord = {}, std::size_t dimension = 0
) {

  // Root case: reached the last dimension
  if (dimension == Tree::dimension) {
    tree.data(tree.child(n, coord.to_ulong())) = points;
    return;
  }

  // Split the point collection around the center point on this dimension
  auto split_point = std::partition(
    points.begin(), points.end(),
    [&](const auto& p) -> bool {
      return (get(point_map, p)[int(dimension)] < center[int(dimension)]);
    }
  );

  // Further subdivide the first side of the split
  std::bitset<Tree::dimension> coord_left = coord;
  coord_left[dimension] = false;
  reassign_points(tree, point_map, n, center, {points.begin(), split_point}, coord_left, dimension + 1);

  // Further subdivide the second side of the split
  std::bitset<Tree::dimension> coord_right = coord;
  coord_right[dimension] = true;
  reassign_points(tree, point_map, n, center, {split_point, points.end()}, coord_right, dimension + 1);
}

/*!
  \ingroup PkgOrthtreeTraits

  Traits class for defining an orthtree of points using the class `CGAL::Orthtree`.

  \tparam GeomTraits model of `Kernel`.
  \tparam PointRange must be a model of `Range` whose value type is the key type of `PointMap` and whose iterator type is model of `RandomAccessIterator`
  \tparam PointMap must be a model of `ReadablePropertyMap` whose value type is a point type from `GeomTraits` matching the current dimension
  \tparam dimension the dimension of the ambient Euclidean space.

  \warning The input point set is not copied. It is used directly
  and is rearranged by the `Orthtree`. Altering the point range
  after creating the orthtree will leave it in an invalid state.

  \cgalModels{CollectionPartitioningOrthtreeTraits}
  \sa `CGAL::Octree`
  \sa `CGAL::Quadtree`
  \sa `CGAL::Orthtree_traits_base<GeomTraits, dimension>`
*/
template <
  typename GeomTraits,
  typename PointRange,
  typename PointMap = Identity_property_map<typename std::iterator_traits<typename PointRange::iterator>::value_type>,
  bool hypercubic_nodes = false,
  int dimension = Ambient_dimension<
    typename std::iterator_traits<typename PointRange::iterator>::value_type,
    GeomTraits
  >::value
>
struct Orthtree_traits_point : public Orthtree_traits_base<GeomTraits, dimension> {
public:
  /// \name Types
  /// @{
  using Node_data = boost::iterator_range<typename PointRange::iterator>;
  /// @}

  using Base = Orthtree_traits_base<GeomTraits, dimension>;
  using Self = Orthtree_traits_point<GeomTraits, PointRange, PointMap, hypercubic_nodes, dimension>;
  using Tree = Orthtree<Self>;

  using Node_index = typename Base::Node_index;
  using Node_data_element = typename std::iterator_traits<typename PointRange::iterator>::value_type;

  Orthtree_traits_point(
    PointRange& points,
    PointMap point_map = PointMap()
  ) : m_points(points), m_point_map(point_map) {}

  auto construct_root_node_bbox_object() const {
    return [&]() -> typename Self::Bbox_d {

      std::array<typename Self::FT, Self::dimension> bbox_min, bbox_max;
      Orthtrees::internal::Cartesian_ranges<Self> cartesian_range;

      // init bbox with first values found
      {
        const typename Self::Point_d& point = get(m_point_map, *(m_points.begin()));
        std::size_t i = 0;
        for (const typename Self::FT& x: cartesian_range(point)) {
          bbox_min[i] = x;
          bbox_max[i] = x;
          ++i;
        }
      }
      // Expand bbox to contain all points
      for (const auto& p: m_points) {
        const typename Self::Point_d& point = get(m_point_map, p);
        std::size_t i = 0;
        for (const typename Self::FT& x: cartesian_range(point)) {
          bbox_min[i] = (std::min)(x, bbox_min[i]);
          bbox_max[i] = (std::max)(x, bbox_max[i]);
          ++i;
        }
      }

#if !defined(_MSC_VER) || _MSC_VER > 1920
      if constexpr (hypercubic_nodes) {
#else
      if (hypercubic_nodes) {
#endif
        std::array<typename Self::FT, Self::dimension> center;
        typename Self::FT max_side = 0;
        for (int i = 0; i < Self::dimension; i++) {
          typename Self::FT side = bbox_max[i] - bbox_min[i];
          max_side = (std::max<typename Self::FT>)(max_side, side);
          center[i] = (bbox_min[i] + bbox_max[i]) * 0.5f;
        }
        max_side *= 0.5f;
        for (int i = 0; i < Self::dimension; i++) {
          bbox_min[i] = center[i] - max_side;
          bbox_max[i] = center[i] + max_side;
        }
      }

      return {std::apply(Self::construct_point_d_object(), bbox_min),
              std::apply(Self::construct_point_d_object(), bbox_max)};
    };
  }

  auto construct_root_node_contents_object() const {
    return [&]() -> typename Self::Node_data {
      return {m_points.begin(), m_points.end()};
    };
  }

  auto distribute_node_contents_object() const {
    return [&](Node_index n, Tree& tree, const typename Self::Point_d& center) {
      CGAL_precondition(!tree.is_leaf(n));
      reassign_points(tree, m_point_map, n, center, tree.data(n));
    };
  }

  auto construct_sphere_d_object() const {
    return [](const typename Self::Point_d& center, const typename Self::FT& squared_radius) -> typename Self::Sphere_d {
      return typename Self::Sphere_d(center, squared_radius);
      };
  }

  auto construct_center_d_object() const {
    return [](const typename Self::Sphere_d& sphere) -> typename Self::Point_d {
      return sphere.center();
      };
  }

  auto compute_squared_radius_d_object() const {
    return [](const typename Self::Sphere_d& sphere) -> typename Self::FT {
      return sphere.squared_radius();
      };
  }

  auto squared_distance_of_element_object() const {
    return [&](const Node_data_element& index, const typename Self::Point_d& point) -> typename Self::FT {
      return CGAL::squared_distance(get(m_point_map, index), point);
      };
  }

  PointRange& m_points;
  PointMap m_point_map;
};

}


#endif //ORTHTREE_TESTS_ORTHTREE_TRAITS_POINT_H
