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

#ifndef ORTHTREE_TESTS_ORTHTREE_TRAITS_POINT_D_H
#define ORTHTREE_TESTS_ORTHTREE_TRAITS_POINT_D_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Dimension.h>
#include <CGAL/Orthtree/Cartesian_ranges.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Kernel_d/Sphere_d.h>

#include <CGAL/Orthtree_traits_d_base.h>

namespace CGAL {

// todo: should this go in its own header & namespace?
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

  The class `Orthtree_traits_point_d` can be used as a template parameter of
  the `Orthtree` class.

  \tparam GeomTraits model of `Kernel`.
  \tparam DimensionTag specialization of `CGAL::Dimension_tag`.
  \tparam PointSet must be a model of range whose value type is the key type of `Point_map`
  \tparam PointMap must be a model of `ReadablePropertyMap` whose value type is `GeomTraits::Traits::Point_d`

  \cgalModels `OrthtreeTraits`
  \sa `CGAL::Orthtree`
  \sa `CGAL::Orthtree_traits_2`
  \sa `CGAL::Orthtree_traits_3`
*/
template <
  typename GeomTraits,
  typename DimensionTag,
  typename PointSet,
  typename PointMap = Identity_property_map<typename GeomTraits::Point_d>
>
struct Orthtree_traits_point_d : public Orthtree_traits_d_base<GeomTraits, DimensionTag> {
public:

  /// \name Types
  /// @{

  using Self = Orthtree_traits_point_d<GeomTraits, DimensionTag, PointSet, PointMap>;

  using Node_data = boost::iterator_range<typename PointSet::iterator>;
  using Node_data_element = typename std::iterator_traits<typename PointSet::iterator>::value_type;

#ifdef DOXYGEN_RUNNING
  /*!
    Functor with an operator to construct a `Point_d` from an `Array` object.
  */
  typedef unspecified_type Construct_point_d_from_array;
#else

  struct Construct_point_d_from_array {
    typename Self::Point_d operator()(const typename Self::Array& array) const {
      return typename Self::Point_d(array.begin(), array.end());
    }
  };

#endif

#ifdef DOXYGEN_RUNNING
  /*!
    Functor with an operator to construct a `Bbox_d` from two `Array` objects (coordinates of minimum and maximum points).
  */
  typedef unspecified_type Construct_bbox_d;
#else

  struct Construct_bbox_d {
    typename Self::Bbox_d operator()(const typename Self::Array& min, const typename Self::Array& max) const {
      return typename Self::Bbox_d(typename Self::Point_d(min.begin(), min.end()), typename Self::Point_d(max.begin(), max.end()));
    }
  };

#endif

  /// @}

  Orthtree_traits_point_d(
    PointSet& point_set,
    PointMap point_map = PointMap()
  ) : m_point_set(point_set), m_point_map(point_map) {}

  /// \name Operations
  /// @{

  /*!
    Function used to construct an object of type `Construct_point_d_from_array`.
  */
  Construct_point_d_from_array construct_point_d_from_array_object() const { return Construct_point_d_from_array(); }

  /*!
    Function used to construct an object of type `Construct_bbox_d`.
  */
  Construct_bbox_d construct_bbox_d_object() const { return Construct_bbox_d(); }

  std::pair<typename Self::Array, typename Self::Array> root_node_bbox() const {

    typename Self::Array bbox_min;
    typename Self::Array bbox_max;
    Orthtrees::internal::Cartesian_ranges<Self> cartesian_range;

    // init bbox with first values found
    {
      const auto& point = get(m_point_map, *(m_point_set.begin()));
      std::size_t i = 0;
      for (const typename Self::FT& x: cartesian_range(point)) {
        bbox_min[i] = x;
        bbox_max[i] = x;
        ++i;
      }
    }
    // Expand bbox to contain all points
    for (const auto& p: m_point_set) {
      const auto& point = get(m_point_map, p);
      std::size_t i = 0;
      for (const typename Self::FT& x: cartesian_range(point)) {
        bbox_min[i] = (std::min)(x, bbox_min[i]);
        bbox_max[i] = (std::max)(x, bbox_max[i]);
        ++i;
      }
    }

    return {bbox_min, bbox_max};
  }

  Node_data root_node_contents() const { return {m_point_set.begin(), m_point_set.end()}; }

  template <typename Node_index, typename Tree>
  void distribute_node_contents(Node_index n, Tree& tree, const typename Self::Point_d& center) {
    CGAL_precondition(!tree.is_leaf(n));
    reassign_points(tree, m_point_map, n, center, tree.data(n));
  }

  typename Self::Point_d get_element(const Node_data_element& index) const {
    return get(m_point_map, index);
  }

  /// @}

private:

  PointSet& m_point_set;
  PointMap m_point_map;

};

}

#endif //ORTHTREE_TESTS_ORTHTREE_TRAITS_POINT_D_H
