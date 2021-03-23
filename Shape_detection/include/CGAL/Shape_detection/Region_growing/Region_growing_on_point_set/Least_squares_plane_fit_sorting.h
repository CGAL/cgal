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

// STL includes.
#include <vector>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>
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
    a model of `LValuePropertyMap` whose key type is the value type of the input
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

    #ifdef DOXYGEN_RUNNING
      /*!
        a model of `ReadablePropertyMap` whose key and value type is `std::size_t`.
        This map provides an access to the ordered indices of input points.
      */
      typedef unspecified_type Seed_map;
    #endif

    /// @}

  private:
    using ITraits = Exact_predicates_inexact_constructions_kernel;
    using IFT = typename ITraits::FT;
    using IPoint_3 = typename ITraits::Point_3;
    using IPlane_3 = typename ITraits::Plane_3;
    using IConverter = Cartesian_converter<Traits, ITraits>;
    using Compare_scores = internal::Compare_scores<IFT>;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

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
    Least_squares_plane_fit_sorting(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      const PointMap point_map = PointMap()) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_point_map(point_map),
    m_iconverter() {

      CGAL_precondition(input_range.size() > 0);
      m_order.resize(m_input_range.size());
      std::iota(m_order.begin(), m_order.end(), 0);
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
    std::vector<std::size_t> m_order;
    std::vector<IFT> m_scores;
    const IConverter m_iconverter;

    void compute_scores() {

      std::vector<IPoint_3> points;
      std::vector<std::size_t> neighbors;
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {

        neighbors.clear();
        m_neighbor_query(i, neighbors);
        neighbors.push_back(i);

        points.clear();
        for (const std::size_t point_index : neighbors) {
          CGAL_precondition(point_index < m_input_range.size());
          const auto& key = *(m_input_range.begin() + point_index);
          const auto& point = get(m_point_map, key);
          points.push_back(m_iconverter(point));
        }
        CGAL_postcondition(points.size() == neighbors.size());

        IPlane_3 fitted_plane;
        IPoint_3 fitted_centroid;
        m_scores[i] = CGAL::linear_least_squares_fitting_3(
          points.begin(), points.end(),
          fitted_plane, fitted_centroid,
          CGAL::Dimension_tag<0>(), ITraits(),
          CGAL::Eigen_diagonalize_traits<IFT, 3>());
      }
    }
  };

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_SORTING_H
