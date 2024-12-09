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

#ifndef CGAL_SHAPE_REGULARIZATION_COLLINEAR_GROUPS_2_H
#define CGAL_SHAPE_REGULARIZATION_COLLINEAR_GROUPS_2_H

#include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/internal/Parallel_groups_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap>
  class Collinear_groups_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;

    using FT = typename Traits::FT;
    using Line_2 = typename Traits::Line_2;
    using Indices = std::vector<std::size_t>;

    using PGroups_2 = Parallel_groups_2<Traits, Input_range, Segment_map>;

    template<typename NamedParameters>
    Collinear_groups_2(
      const InputRange& input_range,
      const NamedParameters& np,
      const SegmentMap segment_map,
      const GeomTraits&) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_grouping(
      input_range, np, segment_map, GeomTraits()) {

      const FT max_offset = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_offset), FT(1) / FT(5));
      const bool preserve_order = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::preserve_order), false);
      CGAL_precondition(max_offset >= FT(0));
      m_max_offset = max_offset;
      make_collinear_groups(preserve_order);
    }

    template<typename OutputIterator>
    OutputIterator groups(OutputIterator groups) const {
      for (const auto& collinear_group : m_collinear_groups) {
        const auto& group = collinear_group;
        *(groups++) = group;
      }
      return groups;
    }

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;
    const PGroups_2 m_grouping;

    FT m_max_offset;
    std::vector<Indices> m_collinear_groups;

    void make_collinear_groups(const bool preserve_order) {

      std::vector<Indices> parallel_groups;
      m_grouping.groups(
        std::back_inserter(parallel_groups));
      m_collinear_groups.reserve(parallel_groups.size());

      Indices collinear_group;
      std::vector<bool> states;

      const FT sq_max_dist = m_max_offset * m_max_offset;
      for (const auto& parallel_group : parallel_groups) {
        CGAL_assertion(parallel_group.size() > 0);

        states.clear();
        states.resize(parallel_group.size(), false);
        handle_parallel_group(
          preserve_order,
          parallel_group, sq_max_dist,
          states, collinear_group);
      }
      CGAL_assertion(
        m_collinear_groups.size() >= parallel_groups.size());
    }

    void handle_parallel_group(
      const bool preserve_order,
      const Indices& parallel_group,
      const FT sq_max_dist,
      std::vector<bool>& states,
      Indices& collinear_group) {

      for (std::size_t i = 0; i < parallel_group.size(); ++i) {
        if (states[i]) continue;

        const std::size_t si_index = parallel_group[i];
        const auto& si = get(m_segment_map,
          *(m_input_range.begin() + si_index));

        states[i] = true;
        collinear_group.clear();
        collinear_group.push_back(si_index);

        const Line_2 line = Line_2(si.source(), si.target());
        traverse_group(
          preserve_order, i, line, parallel_group, sq_max_dist,
          states, collinear_group);
        m_collinear_groups.push_back(collinear_group);
      }
    }

    void traverse_group(
      const bool preserve_order,
      const std::size_t i,
      const Line_2& line,
      const Indices& parallel_group,
      const FT sq_max_dist,
      std::vector<bool>& states,
      Indices& collinear_group) const {

      for (std::size_t j = i + 1; j < parallel_group.size(); ++j) {
        if (states[j]) continue;

        const std::size_t sj_index = parallel_group[j];
        const auto& sj = get(m_segment_map,
          *(m_input_range.begin() + sj_index));

        const auto p = CGAL::midpoint(sj.source(), sj.target());
        const auto q = line.projection(p);

        const FT sq_dist = CGAL::squared_distance(p, q);
        if (sq_dist <= sq_max_dist) {
          states[j] = true;
          collinear_group.push_back(sj_index);
        } else if (preserve_order) return;
      }
    }
  };

} // namespace internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_COLLINEAR_GROUPS_2_H
