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
// Author(s)     : Dmitry Anisimov, Gennadii Sytov
//

#ifndef CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2_H
#define CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2_H

#include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap>
  class Parallel_groups_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;

    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;
    using Indices = std::vector<std::size_t>;

    template<typename NamedParameters>
    Parallel_groups_2(
      const InputRange& input_range,
      const NamedParameters& np,
      const SegmentMap segment_map,
      const GeomTraits&) :
    m_input_range(input_range),
    m_segment_map(segment_map) {

      const FT max_angle = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_angle), FT(5));
      const bool preserve_order = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::preserve_order), false);
      CGAL_precondition(max_angle >= FT(0) && max_angle <= FT(90));
      m_max_angle = max_angle;
      make_parallel_groups(preserve_order);
    }

    template<typename OutputIterator>
    OutputIterator groups(OutputIterator groups) const {
      for (const auto& parallel_group : m_parallel_groups) {
        const auto& group = parallel_group;
        *(groups++) = group;
      }
      return groups;
    }

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;

    FT m_max_angle;
    std::vector<Indices> m_parallel_groups;

    void make_parallel_groups(const bool preserve_order) {

      m_parallel_groups.clear();
      std::vector<bool> states(m_input_range.size(), false);
      Indices parallel_group;

      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        if (states[i]) continue;
        const auto& si = get(
          m_segment_map, *(m_input_range.begin() + i));

        states[i] = true;
        parallel_group.clear();
        parallel_group.push_back(i);

        traverse_group(preserve_order, i, si, states, parallel_group);
        m_parallel_groups.push_back(parallel_group);
      }
    }

    void traverse_group(
      const bool preserve_order,
      const std::size_t i,
      const Segment_2& si,
      std::vector<bool>& states,
      Indices& parallel_group) const {

      for (std::size_t j = i + 1; j < m_input_range.size(); ++j) {
        if (states[j]) continue;
        const auto& sj = get(
          m_segment_map, *(m_input_range.begin() + j));

        const FT angle_2 = internal::angle_2(si, sj);
        if (angle_2 <= m_max_angle) {

          states[j] = true;
          parallel_group.push_back(j);
        } else if (preserve_order) return;
      }
    }
  };

} // namespace internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2_H
