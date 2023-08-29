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

#include <CGAL/Orthtree_traits_point_d.h>
#include <CGAL/Orthtree_traits_2_base.h>

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
struct Orthtree_traits_point_2 : public Orthtree_traits_2_base<GeomTraits> {
public:

  /// \name Types
  /// @{

  using Self = Orthtree_traits_point_2<GeomTraits, PointSet, PointMap>;
  using Tree = Orthtree<Self>;

  // todo: looking for better names
  using Node_data = boost::iterator_range<typename PointSet::iterator>;
  using Node_data_element = typename std::iterator_traits<typename PointSet::iterator>::value_type;

#ifdef DOXYGEN_RUNNING
  /*!
    Functor with an operator to construct a `Point_d` from an `Array` object.
  */
  typedef unspecified_type Construct_point_d_from_array;
#else
  struct Construct_point_d_from_array
  {
    typename Self::Point_d operator() (const typename Self::Array& array) const
    {
      return typename Self::Point_d (array[0], array[1]);
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
    typename Self::Bbox_d operator() (const typename Self::Array& min,
                       const typename Self::Array& max) const
    {
      return typename Self::Bbox_d (min[0], min[1], max[0], max[1]);
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

  auto root_node_bbox_object() const {
    return [&]() -> std::pair<typename Self::Array, typename Self::Array> {

      typename Self::Array bbox_min;
      typename Self::Array bbox_max;
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

      return {bbox_min, bbox_max};
    };
  }

  auto root_node_contents_object() const {
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

  auto get_element_object() const {
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


#endif //ORTHTREE_TESTS_ORTHTREE_TRAITS_POINT_2_H
