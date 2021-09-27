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

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
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

    \brief Fuzzy sphere neighbors search in a set of 2D or 3D points.

    This class returns all neighbors of a query point, which fall in a sphere of
    the fixed radius centered at this point.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Point_2` or `Kernel::Point_3`

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
    using Point_type = typename Point_map::value_type;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// @}

  private:
    using Index_to_point_map =
    internal::Item_property_map<Input_range, Point_map>;

    using Search_base = typename std::conditional<
      std::is_same<typename Traits::Point_2, Point_type>::value,
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

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes a Kd-tree with input points.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param input_range
      an instance of `InputRange` with 2D or 3D points

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{sphere_radius}
          \cgalParamDescription{the fixed radius of the fuzzy sphere used for
          searching neighbors of a query point}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{1}
        \cgalParamNEnd
        \cgalParamNBegin{point_map}
          \cgalParamDescription{an instance of `PointMap` that maps an item from `input_range`
          to `Kernel::Point_2` or to `Kernel::Point_3`}
          \cgalParamDefault{`PointMap()`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() > 0`
      \pre `sphere_radius > 0`
    */
    template<typename NamedParameters>
    Sphere_neighbor_query(
      const InputRange& input_range,
      const NamedParameters& np) :
    m_input_range(input_range),
    m_point_map(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::point_map), PointMap())),
    m_index_to_point_map(m_input_range, m_point_map),
    m_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(m_input_range.size()),
      Splitter(),
      Search_traits(m_index_to_point_map)) {

      CGAL_precondition(input_range.size() > 0);
      m_sphere_radius = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::sphere_radius), FT(1));
      CGAL_precondition(m_sphere_radius > FT(0));
      m_tree.build();
    }

    #if !defined(CGAL_NO_DEPRECATED_CODE) || defined(DOXYGEN_RUNNING)

    /*!
      \brief initializes a Kd-tree with input points.

      \deprecated This constructor is deprecated since the version 5.4 of \cgal.

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
    CGAL_DEPRECATED_MSG("This constructor is deprecated since the version 5.4 of CGAL!")
    Sphere_neighbor_query(
      const InputRange& input_range,
      const FT sphere_radius = FT(1),
      const PointMap point_map = PointMap()) :
    Sphere_neighbor_query(
      input_range, CGAL::parameters::
    sphere_radius(sphere_radius).
    point_map(point_map))
    { }

    #endif // CGAL_NO_DEPRECATED_CODE

    /// \cond SKIP_IN_MANUAL
    Sphere_neighbor_query(
      const InputRange& input_range) :
    Sphere_neighbor_query(
      input_range, CGAL::parameters::all_default())
    { }
    /// \endcond

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `NeighborQuery::operator()()`.

      This operator finds indices of all points, which fall into a sphere
      of the fixed radius `sphere_radius` centered at the query point with
      the index `query_index`. These neighbors are returned in `neighbors`.

      \param query_index
      index of the query point

      \param neighbors
      indices of points, which are neighbors of the query point

      \pre `query_index < input_range.size()`
    */
    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      neighbors.clear();
      CGAL_precondition(query_index < m_input_range.size());
      const std::size_t sphere_center = query_index;
      const Fuzzy_sphere sphere(
        sphere_center, m_sphere_radius, FT(0), m_tree.traits());
      m_tree.search(std::back_inserter(neighbors), sphere);
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    // A property map that can be used to access points when
    // iterating over the indices counting items in the input range.
    const Index_to_point_map& index_to_point_map() const {
      return m_index_to_point_map;
    }
    /// \endcond

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const Index_to_point_map m_index_to_point_map;

    FT m_sphere_radius;
    Tree m_tree;
  };

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_SPHERE_NEIGHBOR_QUERY_H
