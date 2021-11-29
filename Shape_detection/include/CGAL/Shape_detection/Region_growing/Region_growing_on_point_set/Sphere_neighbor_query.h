// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_SPHERE_NEIGHBOR_QUERY_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_SPHERE_NEIGHBOR_QUERY_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <typeinfo>
#include <type_traits>

// Boost includes.
#include <CGAL/boost/iterator/counting_iterator.hpp>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/assertions.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>

namespace CGAL {
namespace Shape_detection {
namespace Point_set {

  /*!
    \ingroup PkgShapeDetectionRGOnPoints

    \brief Fuzzy sphere neighbors search in a set of `Kernel::Point_2`
    or `Kernel::Point_3`.

    This class returns all neighbors of a query point, which fall in a sphere of
    the fixed radius centered at this point.

    \tparam GeomTraits
    must be a model of `Kernel`.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam PointMap
    must be an `LvaluePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Point_2` or `Kernel::Point_3`.

    \cgalModels `NeighborQuery`
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Sphere_neighbor_query {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using Point = typename Point_map::value_type;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// \cond SKIP_IN_MANUAL
    using Index_to_point_map =
    internal::Item_property_map<Input_range, Point_map>;

    using Search_base = typename std::conditional<
      std::is_same<typename Traits::Point_2, Point>::value,
      CGAL::Search_traits_2<Traits>,
      CGAL::Search_traits_3<Traits> >::type;

    using Search_traits =
    CGAL::Search_traits_adapter<std::size_t, Index_to_point_map, Search_base>;

    using Splitter =
    CGAL::Sliding_midpoint<Search_traits>;

    using Fuzzy_sphere
    = CGAL::Fuzzy_sphere<Search_traits>;

    using Tree
    = CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true, CGAL::Tag_true>;
    /// \endcond

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes a Kd-tree with input points.

      \param input_range
      an instance of `InputRange` with 2D or 3D points

      \param sphere_radius
      the fixed radius of the fuzzy sphere used for searching neighbors
      of a query point. %Default is 1.

      \param point_map
      an instance of `PointMap` that maps an item from `input_range`
      to `Kernel::Point_2` or to `Kernel::Point_3`

      \pre `input_range.size() > 0`
      \pre `sphere_radius > 0`
    */
    Sphere_neighbor_query(
      const InputRange& input_range,
      const FT sphere_radius = FT(1),
      const PointMap point_map = PointMap()) :
    m_input_range(input_range),
    m_sphere_radius(sphere_radius),
    m_point_map(point_map),
    m_index_to_point_map(m_input_range, m_point_map),
    m_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(m_input_range.size()),
      Splitter(),
      Search_traits(m_index_to_point_map)) {

      CGAL_precondition(input_range.size() > 0);
      CGAL_precondition(sphere_radius > FT(0));

      m_tree.build();
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `NeighborQuery::operator()()`.

      This operator finds indices of all points, which fall in a sphere
      of the fixed radius `sphere_radius` centered at the query point with
      the index `query_index`. These neighbors are returned in `neighbors`.

      \param query_index
      index of the query point

      \param neighbors
      indices of points, which are neighbors of the query point

      \pre `query_index >= 0 && query_index < input_range.size()`
    */
    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {
      CGAL_precondition(query_index < m_input_range.size());

      const std::size_t sphere_center = query_index;

      const Fuzzy_sphere sphere(
        sphere_center,
        m_sphere_radius,
        FT(0),
        m_tree.traits());

      neighbors.clear();
      m_tree.search(std::back_inserter(neighbors), sphere);
    }

    /// @}

  private:

    // Fields.
    const Input_range& m_input_range;

    const FT m_sphere_radius;

    const Point_map m_point_map;
    const Index_to_point_map m_index_to_point_map;

    Tree m_tree;
  };

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_SPHERE_NEIGHBOR_QUERY_H
