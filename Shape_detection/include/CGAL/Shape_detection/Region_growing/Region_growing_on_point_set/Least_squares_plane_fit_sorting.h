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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_SORTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_SORTING_H

#include <CGAL/license/Shape_detection.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>

namespace CGAL {
namespace Shape_detection {
namespace Point_set {

  /*!
    \ingroup PkgShapeDetectionRGOnPoints

    \brief Sorting of 3D points with respect to the local plane fit quality.

    Indices of 3D input points are sorted with respect to the quality of the
    least squares plane fit applied to the neighboring points of each point.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam NeighborQuery
    a model of `NeighborQuery`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Point_3`
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename NeighborQuery,
  typename PointMap>
  class Least_squares_plane_fit_sorting {

  public:
    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Point_map = PointMap;
    using Seed_map = internal::Seed_property_map;
    /// \endcond

    #ifdef DOXYGEN_NS
      /*!
        a model of `ReadablePropertyMap` whose key and value type is `std::size_t`.
        This map provides an access to the ordered indices of input points.
      */
      typedef unspecified_type Seed_map;
    #endif

    /// @}

  private:
    using FT = typename Traits::FT;
    using Compare_scores = internal::Compare_scores<FT>;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param input_range
      an instance of `InputRange` with 3D points

      \param neighbor_query
      an instance of `NeighborQuery` that is used internally to
      access point's neighbors

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{point_map}
          \cgalParamDescription{an instance of `PointMap` that maps an item from `input_range`
          to `Kernel::Point_3`}
          \cgalParamDefault{`PointMap()`}
        \cgalParamNEnd
        \cgalParamNBegin{geom_traits}
          \cgalParamDescription{an instance of `GeomTraits`}
          \cgalParamDefault{`GeomTraits()`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() > 0`
    */
    template<typename NamedParameters>
    Least_squares_plane_fit_sorting(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      const NamedParameters& np) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_point_map(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::point_map), PointMap())),
    m_traits(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::geom_traits), GeomTraits())) {

      CGAL_precondition(input_range.size() > 0);
      m_order.resize(m_input_range.size());
      std::iota(m_order.begin(), m_order.end(), 0);
      m_scores.resize(m_input_range.size());
    }

    /*!
      \brief initializes all internal data structures.

      \deprecated This constructor is deprecated since the version 5.4 of \cgal.

      \param input_range
      an instance of `InputRange` with 3D points

      \param neighbor_query
      an instance of `NeighborQuery` that is used internally to
      access point's neighbors

      \param point_map
      an instance of `PointMap` that maps an item from `input_range`
      to `Kernel::Point_3`

      \pre `input_range.size() > 0`
    */
    CGAL_DEPRECATED_MSG("This constructor is deprecated since the version 5.4 of CGAL!")
    Least_squares_plane_fit_sorting(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      const PointMap point_map = PointMap()) :
    Least_squares_plane_fit_sorting(
      input_range, neighbor_query, CGAL::parameters::
    point_map(point_map))
    { }

    /// @}

    /// \name Sorting
    /// @{

    /*!
      \brief sorts indices of input points.
    */
    void sort() {

      compute_scores();
      CGAL_precondition(m_scores.size() > 0);
      Compare_scores cmp(m_scores);
      std::sort(m_order.begin(), m_order.end(), cmp);
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief returns an instance of `Seed_map` to access the ordered indices
      of input points.
    */
    Seed_map seed_map() {
      return Seed_map(m_order);
    }

    /// @}

  private:
    const Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    const Point_map m_point_map;
    const Traits m_traits;
    std::vector<std::size_t> m_order;
    std::vector<FT> m_scores;

    void compute_scores() {

      std::vector<std::size_t> neighbors;
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        neighbors.clear();
        m_neighbor_query(i, neighbors);
        neighbors.push_back(i);
        m_scores[i] = internal::create_plane(
          m_input_range, m_point_map, neighbors, m_traits).second;
      }
    }
  };

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_SORTING_H
