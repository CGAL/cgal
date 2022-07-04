// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_LEAST_SQUARES_LINE_FIT_SORTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_LEAST_SQUARES_LINE_FIT_SORTING_H

#include <CGAL/license/Shape_detection.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>

namespace CGAL {
namespace Shape_detection {
namespace Polyline {

  /*!
    \ingroup PkgShapeDetectionRGOnPolyline

    \brief Sorting of polyline vertices with respect to the local line fit quality.

    Indices of input vertices are sorted with respect to the quality of the
    least squares line fit applied to the incident vertices of each vertex.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam NeighborQuery
    a model of `NeighborQuery`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Point_2` or `Kernel::Point_3`
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename NeighborQuery,
  typename PointMap>
  class Least_squares_line_fit_sorting {

  public:
    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Point_map = PointMap;
    using Point_type = typename Point_map::value_type;
    /// \endcond

    /// Item type.
    using Item = typename InputRange::const_iterator;

    /// Seed range.
    using Seed_range = std::vector<Item>;

    /// @}

  private:
    using FT = typename Traits::FT;
    using Polyline_traits = typename std::conditional<
      std::is_same<typename Traits::Point_2, Point_type>::value,
      internal::Region_growing_traits_2<Traits>,
      internal::Region_growing_traits_3<Traits> >::type;
    using Compare_scores = internal::Compare_scores<FT>;
    using Dereference_pmap = internal::Dereference_property_map_adaptor<Item, PointMap>;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param input_range
      an instance of `InputRange` with polyline vertices

      \param neighbor_query
      an instance of `NeighborQuery` that is used internally to
      access vertex's neighbors

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{point_map}
          \cgalParamDescription{an instance of `PointMap` that maps a vertex from `input_range`
          to `Kernel::Point_2` or `Kernel::Point_3`}
          \cgalParamDefault{`PointMap()`}
        \cgalParamNEnd
        \cgalParamNBegin{geom_traits}
          \cgalParamDescription{an instance of `GeomTraits`}
          \cgalParamDefault{`GeomTraits()`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() > 0`
    */
    template<typename NamedParameters = parameters::Default_named_parameters>
    Least_squares_line_fit_sorting(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      const NamedParameters& np = parameters::default_values()) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_point_map(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::point_map), PointMap())),
    m_deref_pmap(m_point_map),
    m_traits(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::geom_traits), GeomTraits())),
    m_polyline_traits(m_traits) {

      CGAL_precondition(input_range.size() > 0);

      m_ordered.resize(m_input_range.size());

      std::size_t index = 0;
      for (auto it = m_input_range.begin(); it != m_input_range.end(); it++)
        m_ordered[index++] = it;

      m_scores.resize(m_input_range.size());
    }

    /// @}

    /// \name Sorting
    /// @{

    /*!
      \brief sorts `Items` of input vertices.
    */
    void sort() {
      compute_scores();
      CGAL_precondition(m_scores.size() > 0);
      Compare_scores cmp(m_scores);

      std::vector<std::size_t> order(m_input_range.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), cmp);

      std::vector<Item> tmp(m_input_range.size());
      for (std::size_t i = 0; i < m_input_range.size(); i++)
        tmp[i] = m_ordered[order[i]];

      m_ordered.swap(tmp);
    }
    /// @}

    /// \name Access
    /// @{

    /*!
      \brief returns an instance of `Seed_range` to access the ordered `Items`
      of input points.
    */
      const Seed_range& ordered() {
      return m_ordered;
    }
    /// @}

  private:
    const Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    const Point_map m_point_map;
    Dereference_pmap m_deref_pmap;
    const Traits m_traits;
    const Polyline_traits m_polyline_traits;
    Seed_range m_ordered;
    std::vector<FT> m_scores;

    void compute_scores() {

      std::vector<Item> neighbors;
      std::size_t idx = 0;
      for (auto it = m_input_range.begin(); it != m_input_range.end(); it++) {
        neighbors.clear();
        m_neighbor_query(it, neighbors);
        neighbors.push_back(it);
        m_scores[idx++] = m_polyline_traits.create_line(
          neighbors, m_point_map).second;
      }
    }
  };

} // namespace Polyline
} // namespace Shape_detection
} // namespace CGAL

#endif // #define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_LEAST_SQUARES_LINE_FIT_SORTING_H
