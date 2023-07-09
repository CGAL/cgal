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
// Author(s)     : Dmitry Anisimov, Simon Giraudot
//

#ifndef CGAL_SHAPE_REGULARIZATION_LONGEST_PRINCIPAL_DIRECTION_2_H
#define CGAL_SHAPE_REGULARIZATION_LONGEST_PRINCIPAL_DIRECTION_2_H

#include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Contour_base_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace Contours {

  /*!
    \ingroup PkgShapeRegularizationRefContours

    \brief Estimates the longest principal direction of the contour.

    This algorithm finds the longest contour edge and sets its direction as the
    principal direction of the contour.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `GeomTraits::Point_2`. The default is
    `CGAL::Identity_property_map<typename GeomTraits::Point_2>`.

    \cgalModels{ContourDirections}
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Longest_direction_2 {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;
    using Direction_2 = typename Traits::Direction_2;

    using FT_pair = std::pair<FT, FT>;
    using Base = internal::Contour_base_2<Traits>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param input_range
      a const range of ordered 2D points, which form a contour

      \param is_closed
      indicates whether the contour is closed or open

      \param np
      an optional sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below; this parameter can be omitted,
      the default values are then used

      \cgalNamedParamsBegin
        \cgalParamNBegin{point_map}
          \cgalParamDescription{an instance of `PointMap` that maps an item from `input_range`
          to `GeomTraits::Point_2`}
          \cgalParamDefault{`PointMap()`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre input_range.size() >= 3 for closed contours
      \pre input_range.size() >= 2 for open contours
    */
    template<typename NamedParameters = parameters::Default_named_parameters>
    Longest_direction_2(
      const InputRange& input_range,
      const bool is_closed,
      const NamedParameters& np = parameters::default_values()) :
    m_input_range(input_range),
    m_point_map(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::point_map), PointMap())) {

      CGAL_precondition(m_input_range.size() >= 2);

      if (is_closed) {
        estimate_closed(m_bounds, m_directions, m_assigned);
      } else {
        estimate_open(m_bounds, m_directions, m_assigned);
      }

      if (verbose()) {
        std::cout << "* assigned directions: ";
        for (std::size_t direction_index : m_assigned) {
          std::cout << direction_index << " ";
        }
        std::cout << std::endl;
      }
    }

    /// @}

    /// \name Directions
    /// @{

    /*!
      \brief orients a given `segment` with the index `query_index` towards the
      longest principal direction.

      \param query_index
      an index of the contour vertex that emits the contour edge being `segment`

      \param segment
      a segment to be rotated

      \pre query_index >= 0 && query_index < input_range.size() for closed contours
      \pre query_index >= 0 && query_index < input_range.size() - 1 for open contours
    */
    void orient(
      const std::size_t query_index,
      Segment_2& segment) const {

      m_base.apply_rotation_to_segment(
        m_bounds, m_directions, m_assigned,
        query_index, segment);
    }

    /// @}

    /// \name Miscellaneous
    /// @{

    /*!
      \brief returns the number of principal directions of the contour.

      The returned number is always one.
    */
    std::size_t number_of_directions() const {
      return m_directions.size();
    }

    /// @}

    // EXTRA METHODS TO TEST THE CLASS!
    /// \cond SKIP_IN_MANUAL
    const std::vector<Direction_2>& get_directions() const {
      return m_directions;
    }
    /// \endcond

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const Base m_base;

    std::vector<FT_pair> m_bounds;
    std::vector<Direction_2> m_directions;
    std::vector<std::size_t> m_assigned;

    bool verbose() const {
      return m_base.verbose();
    }

    void estimate_closed(
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      bounds.clear(); bounds.resize(1);
      bounds[0] = std::make_pair(FT(45), FT(45));

      directions.clear(); directions.resize(1);
      directions[0] = compute_longest_direction_closed();

      // 0 is the index of the direction in the `directions`.
      assigned.clear();
      assigned.resize(m_input_range.size(), 0);
    }

    Direction_2 compute_longest_direction_closed() const {

      FT max_length = -FT(1);
      std::size_t index = std::size_t(-1);

      CGAL_precondition(m_input_range.size() >= 3);
      const std::size_t n = m_input_range.size();
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ip = (i + 1) % n;

        const auto& source = get(m_point_map, *(m_input_range.begin() + i));
        const auto& target = get(m_point_map, *(m_input_range.begin() + ip));

        const FT sq_length = CGAL::squared_distance(source, target);
        if (sq_length > max_length) {
          index = i; max_length = sq_length;
        }
      }

      const std::size_t i = index;
      const std::size_t ip = (i + 1) % n;
      const auto& source = get(m_point_map, *(m_input_range.begin() + i));
      const auto& target = get(m_point_map, *(m_input_range.begin() + ip));

      const Segment_2 segment = Segment_2(source, target);
      auto v = segment.to_vector();
      return internal::direction_2(v);
    }

    void estimate_open(
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      bounds.clear(); bounds.resize(1);
      bounds[0] = std::make_pair(FT(45), FT(45));

      directions.clear(); directions.resize(1);
      directions[0] = compute_longest_direction_open();

      // 0 is the index of the direction in the `directions`.
      assigned.clear();
      assigned.resize(m_input_range.size() - 1, 0);
    }

    Direction_2 compute_longest_direction_open() const {

      FT max_length = -FT(1);
      std::size_t index = std::size_t(-1);

      CGAL_assertion(m_input_range.size() >= 2);
      const std::size_t n = m_input_range.size();
      for (std::size_t i = 0; i < n - 1; ++i) {
        const std::size_t ip = i + 1;

        const auto& source = get(m_point_map, *(m_input_range.begin() + i));
        const auto& target = get(m_point_map, *(m_input_range.begin() + ip));

        const FT sq_length = CGAL::squared_distance(source, target);
        if (sq_length > max_length) {
          index = i; max_length = sq_length;
        }
      }

      const std::size_t i = index;
      const std::size_t ip = i + 1;
      const auto& source = get(m_point_map, *(m_input_range.begin() + i));
      const auto& target = get(m_point_map, *(m_input_range.begin() + ip));

      const Segment_2 segment = Segment_2(source, target);
      auto v = segment.to_vector();
      return internal::direction_2(v);
    }
  };

} // namespace Contours
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_LONGEST_PRINCIPAL_DIRECTION_2_H
