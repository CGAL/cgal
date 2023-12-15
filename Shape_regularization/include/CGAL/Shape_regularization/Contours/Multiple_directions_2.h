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

#ifndef CGAL_SHAPE_REGULARIZATION_MULTIPLE_PRINCIPAL_DIRECTIONS_2_H
#define CGAL_SHAPE_REGULARIZATION_MULTIPLE_PRINCIPAL_DIRECTIONS_2_H

#include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Contour_base_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace Contours {

  /*!
    \ingroup PkgShapeRegularizationRefContours

    \brief Estimates possibly multiple principal directions of the contour
    based on the user-specified minimum length and maximum angle bounds.

    This algorithm finds the best-fit edges of the contour with respect to the
    user-specified parameters and sets their directions as the principal directions
    of the contour.

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
  class Multiple_directions_2 {

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
    using Segment_wrapper_2 = typename Base::Segment_wrapper_2;
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
        \cgalParamNBegin{maximum_angle}
          \cgalParamDescription{maximum allowed angle deviation in degrees between a contour edge
            and a principal direction such that they are considered to be parallel or orthogonal}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{10 degrees}
        \cgalParamNEnd
        \cgalParamNBegin{minimum_length}
          \cgalParamDescription{minimum acceptable length of a contour edge whose direction can be taken
            as a principal direction}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{3 unit lengths}
        \cgalParamNEnd
        \cgalParamNBegin{adjust_directions}
          \cgalParamDescription{indicates whether the found directions should be better
            adjusted to the original shape or not}
          \cgalParamType{boolean}
          \cgalParamDefault{true}
        \cgalParamNEnd
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
    Multiple_directions_2(
      const InputRange& input_range,
      const bool is_closed,
      const NamedParameters& np = parameters::default_values()) :
    m_input_range(input_range),
    m_point_map(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::point_map), PointMap())) {

      CGAL_precondition(input_range.size() >= 2);
      m_max_angle_2 = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_angle), FT(10));
      m_min_length_2 = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::minimum_length), FT(3));
      m_adjust_directions = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::adjust_directions), true);

      CGAL_precondition(
        m_max_angle_2 >= FT(0) &&
        m_max_angle_2 <= FT(90));
      CGAL_precondition(
        m_min_length_2 >= FT(0));

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
      best-fit found principal direction.

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

      The returned number is always greater or equal than one.
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

    FT m_max_angle_2;
    FT m_min_length_2;
    bool m_adjust_directions;

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

      std::vector<Segment_wrapper_2> wraps;
      m_base.initialize_closed(
        m_input_range, m_point_map, wraps);

      if (wraps.size() < 3) {
        set_longest_direction(
          wraps, bounds, directions, assigned);
        return;
      }

      set_valid_directions(wraps);
      estimate_initial_directions(
        wraps, bounds, directions, assigned);

      if (directions.size() <= 1) {
        set_longest_direction(
          wraps, bounds, directions, assigned);
      } else {
        m_base.unify_along_contours_closed(wraps, assigned);
        m_base.correct_directions_closed(wraps, assigned);
        if (m_adjust_directions) {
          m_base.readjust_directions(wraps, assigned, directions);
        }
      }
    }

    void estimate_open(
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      std::vector<Segment_wrapper_2> wraps;
      m_base.initialize_open(
        m_input_range, m_point_map, wraps);

      if (wraps.size() < 3) {
        set_longest_direction(
          wraps, bounds, directions, assigned);
        return;
      }

      set_valid_directions(wraps);
      estimate_initial_directions(
        wraps, bounds, directions, assigned);

      if (directions.size() <= 1) {
        set_longest_direction(
          wraps, bounds, directions, assigned);
      } else {
        m_base.unify_along_contours_open(wraps, assigned);
        m_base.correct_directions_open(wraps, assigned);
        if (m_adjust_directions) {
          m_base.readjust_directions(wraps, assigned, directions);
        }
      }
    }

    void set_valid_directions(
      std::vector<Segment_wrapper_2>& wraps) const {

      for (auto& wrap : wraps) {
        wrap.is_valid_direction =
          is_valid_principal_direction(wrap.segment);
      }
    }

    bool is_valid_principal_direction(
      const Segment_2& segment) const {

      CGAL_assertion(m_min_length_2 >= FT(0));
      const FT threshold = m_min_length_2 * FT(2);
      const FT squared_threshold = threshold * threshold;
      return segment.squared_length() >= squared_threshold;
    }

    void estimate_initial_directions(
      std::vector<Segment_wrapper_2>& wraps,
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      std::vector<std::size_t> longest_to_short;
      m_base.sort_segments_by_length(wraps, longest_to_short);
      CGAL_assertion(longest_to_short.size() == wraps.size());

      bounds.clear(); directions.clear(); assigned.clear();
      assigned.resize(longest_to_short.size(), std::size_t(-1));

      std::size_t group_index = 0;
      std::size_t query_index = std::size_t(-1);
      do {
        query_index = find_next_longest_segment(
          wraps, longest_to_short);
        if (query_index != std::size_t(-1)) {
          set_next_longest_direction(
            wraps, query_index, group_index,
            bounds, directions, assigned);
        }
        ++group_index;
      } while (query_index != std::size_t(-1));
    }

    std::size_t find_next_longest_segment(
      const std::vector<Segment_wrapper_2>& wraps,
      const std::vector<std::size_t>& longest_to_short) const {

      std::size_t longest = std::size_t(-1);
      for (std::size_t i = 0; i < longest_to_short.size(); ++i) {
        const std::size_t wrap_index = longest_to_short[i];
        const auto& wrap = wraps[wrap_index];
        if (is_valid_wrap(wrap)) {
          longest = wrap_index; break;
        }
      }
      return longest;
    }

    bool is_valid_wrap(
      const Segment_wrapper_2& wrap) const {
      return !wrap.is_used && wrap.is_valid_direction;
    }

    void set_next_longest_direction(
      std::vector<Segment_wrapper_2>& wraps,
      const std::size_t query_index,
      const std::size_t group_index,
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      CGAL_assertion(query_index != std::size_t(-1));
      CGAL_assertion(group_index != std::size_t(-1));

      // Set current longest direction.
      auto& longest = wraps[query_index];
      assigned[query_index] = group_index;
      longest.is_used = true;

      for (auto& wrap : wraps) {
        if (wrap.index == query_index) { // skip longest
          continue;
        }

        // Check if another wrap satisfies the conditions.
        if (is_valid_wrap(wrap)) {
          if (does_satisify_angle_conditions(
            longest.segment, wrap.segment)) {

            assigned[wrap.index] = group_index;
            wrap.is_used = true;
          }
        }
      }

      // Set internals.
      directions.push_back(longest.direction);
      bounds.push_back(std::make_pair(FT(45), FT(45)));
    }

    bool does_satisify_angle_conditions(
      const Segment_2& longest,
      const Segment_2& segment) const {

      CGAL_precondition(
        m_max_angle_2 >= FT(0) && m_max_angle_2 <= FT(90));
      const FT bound_min = m_max_angle_2;
      const FT bound_max = FT(90) - bound_min;

      const FT angle_2 = internal::angle_2(longest, segment);
      return (angle_2 <= bound_min) || (angle_2 >= bound_max);
    }

    void set_longest_direction(
      const std::vector<Segment_wrapper_2>& wraps,
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      bounds.clear(); bounds.resize(1);
      bounds[0] = std::make_pair(FT(45), FT(45));

      directions.clear(); directions.resize(1);
      directions[0] = compute_longest_direction(wraps);

      // 0 is the index of the direction in the `directions`.
      assigned.clear();
      assigned.resize(wraps.size(), 0);
    }

    Direction_2 compute_longest_direction(
      const std::vector<Segment_wrapper_2>& wraps) const {

      const std::size_t n = wraps.size();
      CGAL_assertion(n != 0);

      FT max_length = -FT(1);
      std::size_t longest = std::size_t(-1);

      for (std::size_t i = 0; i < n; ++i) {
        const auto& wrap = wraps[i];
        const FT sq_length = wrap.segment.squared_length();
        if (sq_length > max_length) {
          longest = i; max_length = sq_length;
        }
      }

      CGAL_assertion(longest != std::size_t(-1));
      return wraps[longest].direction;
    }
  };

} // namespace Contours
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_MULTIPLE_PRINCIPAL_DIRECTIONS_2_H
