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

#ifndef ORTHTREE_TESTS_ORTHTREE_TRAITS_POINT_2_H
#define ORTHTREE_TESTS_ORTHTREE_TRAITS_POINT_2_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Dimension.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Point_set_2.h>
#include <CGAL/Orthtree/Cartesian_ranges.h>

namespace CGAL {

/*!
  \ingroup PkgOrthtreeTraits

  The class `Orthtree_traits_point_2` can be used as a template parameter of
  the `Orthtree` class.

  \tparam GeomTraits model of `Kernel`.
  \tparam PointSet must be a model of range whose value type is the key type of `PointMap`
  \tparam PointMap must be a model of `ReadablePropertyMap` whose value type is `GeomTraits::Traits::Point_d`

  \cgalModels `OrthtreeTraits`
  \sa `CGAL::Octree`
  \sa `CGAL::Orthtree_traits_2`
  \sa `CGAL::Orthtree_traits_d`
*/
template <
  typename GeomTraits,
  typename PointSet,
  typename PointMap = Identity_property_map<typename GeomTraits::Point_2>
>
struct Orthtree_traits_point_2 {
public:

  /// \name Types
  /// @{

  using Self = Orthtree_traits_point_2<GeomTraits, PointSet, PointMap>;

  using Dimension = Dimension_tag<2>;
  using Bbox_d = Bbox_2;
  using FT = typename GeomTraits::FT;
  using Point_d = typename GeomTraits::Point_2;
  using Sphere_d = typename GeomTraits::Circle_2;
  using Cartesian_const_iterator_d = typename GeomTraits::Cartesian_const_iterator_2;
  using Array = std::array<FT, Dimension::value>; // todo: This should have a more descriptive name

  // todo: looking for better names
  using Node_data = boost::iterator_range<typename PointSet::iterator>;
  using Node_data_element = typename std::iterator_traits<typename PointSet::iterator>::value_type;

  /*!
   * \brief Two directions along each axis in Cartesian space, relative to a node.
   *
   * Directions are mapped to numbers as 2-bit integers.
   *
   * The first bit indicates the axis (0 = x, 1 = y),
   * the second bit indicates the direction along that axis (0 = -, 1 = +).
   *
   * The following diagram may be a useful reference:
   *
   *            3 *
   *              |
   *              |                    y+
   *              |                    *
   *     0 *------+------* 1           |
   *              |                    |
   *              |                    +-----* x+
   *              |
   *              * 2
   *
   * This lookup table may also be helpful:
   *
   * | Direction | bitset | number | Enum  |
   * | --------- | ------ | ------ | ----- |
   * | `-x`      |  00    | 0      | LEFT  |
   * | `+x`      |  01    | 1      | RIGHT |
   * | `-y`      |  10    | 2      | DOWN  |
   * | `+y`      |  11    | 3      | UP    |
   */
  enum Adjacency
  {
    LEFT,
    RIGHT,
    DOWN,
    UP
  };

#ifdef DOXYGEN_RUNNING
  /*!
    Functor with an operator to construct a `Point_d` from an `Array` object.
  */
  typedef unspecified_type Construct_point_d_from_array;
#else
  struct Construct_point_d_from_array
  {
    Point_d operator() (const Array& array) const
    {
      return Point_d (array[0], array[1]);
    }
  };
#endif


#ifdef DOXYGEN_RUNNING
  /*!
    Functor with an operator to construct a `Bbox_d` from two `Array` objects (coordinates of minimum and maximum points).
  */
  typedef unspecified_type Construct_bbox_d;
#else
  struct Construct_bbox_d
  {
    Bbox_d operator() (const Array& min,
                       const Array& max) const
    {
      return Bbox_d (min[0], min[1], max[0], max[1]);
    }
  };
#endif

  /// @}

  Orthtree_traits_point_2(
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

  std::pair<Array, Array> root_node_bbox() const {

    Array bbox_min;
    Array bbox_max;
    Orthtrees::internal::Cartesian_ranges<Self> cartesian_range;

    // init bbox with first values found
    {
      const Point_d& point = get(m_point_map, *(m_point_set.begin()));
      std::size_t i = 0;
      for (const FT& x: cartesian_range(point)) {
        bbox_min[i] = x;
        bbox_max[i] = x;
        ++i;
      }
    }
    // Expand bbox to contain all points
    for (const auto& p: m_point_set) {
      const Point_d& point = get(m_point_map, p);
      std::size_t i = 0;
      for (const FT& x: cartesian_range(point)) {
        bbox_min[i] = (std::min)(x, bbox_min[i]);
        bbox_max[i] = (std::max)(x, bbox_max[i]);
        ++i;
      }
    }

    return {bbox_min, bbox_max};
  }

  Node_data root_node_contents() const { return {m_point_set.begin(), m_point_set.end()}; }

  template <typename Node_index, typename Tree>
  void distribute_node_contents(Node_index n, Tree& tree, const Point_d& center) {
    CGAL_precondition(!tree.is_leaf(n));
    reassign_points(n, tree, center, tree.data(n));
  }

  Point_d get_element(const Node_data_element& index) const {
    return get(m_point_map, index);
  }

  /// @}

private:

  PointSet& m_point_set;
  PointMap m_point_map;

  template <typename Node_index, typename Tree>
  void reassign_points(Node_index n, Tree& tree,
                       const Point_d& center, Node_data points,
                       std::bitset<Dimension::value> coord = {},
                       std::size_t dimension = 0) {

    // Root case: reached the last dimension
    if (dimension == Dimension::value) {
      tree.data(tree.child(n, coord.to_ulong())) = points;
      return;
    }

    // Split the point collection around the center point on this dimension
    auto split_point = std::partition(
      points.begin(), points.end(),
      [&](const auto& p) -> bool {
        // This should be done with cartesian iterator,
        // but it seems complicated to do efficiently
        return (get(m_point_map, p)[int(dimension)] < center[int(dimension)]);
      }
    );

    // Further subdivide the first side of the split
    std::bitset<Dimension::value> coord_left = coord;
    coord_left[dimension] = false;
    reassign_points(n, tree, center, {points.begin(), split_point}, coord_left, dimension + 1);

    // Further subdivide the second side of the split
    std::bitset<Dimension::value> coord_right = coord;
    coord_right[dimension] = true;
    reassign_points(n, tree, center, {split_point, points.end()}, coord_right, dimension + 1);
  }

};

}


#endif //ORTHTREE_TESTS_ORTHTREE_TRAITS_POINT_2_H
