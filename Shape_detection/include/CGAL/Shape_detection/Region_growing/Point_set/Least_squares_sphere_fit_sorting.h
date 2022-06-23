// Copyright (c) 2021 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_SPHERE_FIT_SORTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_SPHERE_FIT_SORTING_H

#include <CGAL/license/Shape_detection.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>

namespace CGAL {
namespace Shape_detection {
namespace Point_set {

  /*!
    \ingroup PkgShapeDetectionRGOnPoints

    \brief Sorting of 3D points with respect to the local sphere fit quality.

    Indices of 3D input points are sorted with respect to the quality of the
    least squares sphere fit applied to the neighboring points of each point.

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
  class Least_squares_sphere_fit_sorting {

  public:
    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Point_map = PointMap;

    using Item = typename InputRange::const_iterator;
    using Region = std::vector<Item>;
    using Seed_range = std::vector<Item>;
    /// \endcond

    #ifdef DOXYGEN_RUNNING
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
    template<typename CGAL_NP_TEMPLATE_PARAMETERS>
    Least_squares_sphere_fit_sorting(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      const CGAL_NP_CLASS& np = parameters::default_values()) :
      m_input_range(input_range),
      m_neighbor_query(neighbor_query),
      m_point_map(Point_set_processing_3_np_helper<InputRange, CGAL_NP_CLASS,PointMap>::get_const_point_map(input_range, np)),
      m_traits(parameters::choose_parameter(parameters::get_parameter(
        np, internal_np::geom_traits), GeomTraits())) {

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
      \brief sorts indices of input points.
    */
    void sort() {

      compute_scores();
      CGAL_postcondition(m_scores.size() > 0);
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
      of input points.
    */
    const Seed_range &ordered() {
      return m_ordered;
    }

    /// @}

  private:
    const Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    const Point_map m_point_map;
    const Traits m_traits;
    Seed_range m_ordered;
    std::vector<FT> m_scores;

    void compute_scores() {

      std::vector<Item> neighbors;
      std::size_t idx = 0;
      for (auto it = m_input_range.begin(); it != m_input_range.end(); it++) {
        neighbors.clear();
        m_neighbor_query(it, neighbors);
        neighbors.push_back(it);

        m_scores[idx] = internal::create_sphere(
          neighbors, m_point_map, m_traits, true).second;
        CGAL_assertion(m_scores[idx] <= FT(0));
        idx++;
      }
    }
  };

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_SPHERE_FIT_SORTING_H
