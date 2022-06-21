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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LEAST_SQUARES_LINE_FIT_SORTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LEAST_SQUARES_LINE_FIT_SORTING_H

#include <CGAL/license/Shape_detection.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>

namespace CGAL {
namespace Shape_detection {
namespace Segment_set {

  /*!
    \ingroup PkgShapeDetectionRGOnSegments

    \brief Sorting of segments with respect to the local line fit quality.

    Indices of input segments are sorted with respect to the quality of the
    least squares line fit applied to the vertices of incident segments of each segment.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam NeighborQuery
    a model of `NeighborQuery`

    \tparam SegmentMap
    a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Segment_2` or `Kernel::Segment_3`
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename NeighborQuery,
  typename SegmentMap>
  class Least_squares_line_fit_sorting {

  public:
    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Segment_map = SegmentMap;
    using Segment_type = typename Segment_map::value_type;

    using Item = typename InputRange::const_iterator;
    using Region = std::vector<Item>;
    using Seed_range = std::vector<Item>;
    /// \endcond

    #ifdef DOXYGEN_NS
      /*!
        a model of `ReadablePropertyMap` whose key and value type is `std::size_t`.
        This map provides an access to the ordered indices of input segments.
      */
      typedef unspecified_type Seed_map;
    #endif

    /// @}

  private:
    using FT = typename Traits::FT;
    using Segment_set_traits = typename std::conditional<
      std::is_same<typename Traits::Segment_2, Segment_type>::value,
      internal::Region_growing_traits_2<Traits>,
      internal::Region_growing_traits_3<Traits> >::type;
    using Compare_scores = internal::Compare_scores<FT>;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param input_range
      an instance of `InputRange` with 2D or 3D segments

      \param neighbor_query
      an instance of `NeighborQuery` that is used internally to
      access vertex's neighbors

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{segment_map}
          \cgalParamDescription{an instance of `SegmentMap` that maps a segment from `input_range`
          to `Kernel::Segment_2` or `Kernel::Segment_3`}
          \cgalParamDefault{`SegmentMap()`}
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
    m_segment_map(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::segment_map), SegmentMap())),
    m_traits(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::geom_traits), GeomTraits())),
    m_segment_set_traits(m_traits) {

      CGAL_precondition(input_range.size() > 0);

      m_ordered.resize(m_input_range.size());

      std::size_t index = 0;
      for (Item& i = m_input_range.begin(); i != m_input_range.end(); i++)
        m_ordered[index++] = i;

      m_scores.resize(m_input_range.size());
    }

    /// @}

    /// \name Sorting
    /// @{

    /*!
      \brief sorts indices of input segments.
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
      \brief returns an instance of `Seed_map` to access the ordered indices
      of input segments.
    */
    const Seed_range &ordered() {
      return m_ordered;
    }

    /// @}

  private:
    const Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    const Segment_map m_segment_map;
    const Traits m_traits;
    const Segment_set_traits m_segment_set_traits;
    Seed_range m_ordered;
    std::vector<FT> m_scores;

    void compute_scores() {

      std::vector<Item> neighbors;
      std::size_t idx = 0;
      for (auto it = m_input_range.begin(); it != m_input_range.end(); it++) {
        neighbors.clear();
        m_neighbor_query(it, neighbors);
        neighbors.push_back(it);

        const auto& segment = get(m_segment_map, *it);
        const auto& source = segment.source();
        const auto& target = segment.target();
        if (source == target)
          m_scores[idx++] = FT(0); // put it at the very back
        else
          m_scores[idx++] =
            m_segment_set_traits.create_line(neighbors, m_segment_map).second;
      }
    }
  };

} // namespace Segment_set
} // namespace Shape_detection
} // namespace CGAL

#endif // #define CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LEAST_SQUARES_LINE_FIT_SORTING_H
