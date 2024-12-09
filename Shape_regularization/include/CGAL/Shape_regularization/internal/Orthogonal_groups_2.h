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

#ifndef CGAL_SHAPE_REGULARIZATION_ORTHOGONAL_GROUPS_2_H
#define CGAL_SHAPE_REGULARIZATION_ORTHOGONAL_GROUPS_2_H

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
  class Orthogonal_groups_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;

    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;
    using Indices = std::vector<std::size_t>;

    using PGroups_2 = Parallel_groups_2<Traits, Input_range, Segment_map>;

    template<typename NamedParameters>
    Orthogonal_groups_2(
      const InputRange& input_range,
      const NamedParameters& np,
      const SegmentMap segment_map,
      const GeomTraits&) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_grouping(
      input_range, np, segment_map, GeomTraits()) {

      const FT max_angle = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_angle), FT(5));
      const bool preserve_order = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::preserve_order), false);
      CGAL_precondition(max_angle >= FT(0) && max_angle <= FT(90));
      m_max_angle = max_angle;
      make_orthogonal_groups(preserve_order);
    }

    template<typename OutputIterator>
    OutputIterator groups(OutputIterator groups) const {
      for (const auto& orthogonal_group : m_orthogonal_groups) {
        const auto& group = orthogonal_group;
        *(groups++) = group;
      }
      return groups;
    }

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;
    const PGroups_2 m_grouping;

    FT m_max_angle;
    std::vector<Indices> m_orthogonal_groups;

    void make_orthogonal_groups(const bool preserve_order) {

      std::vector<Indices> parallel_groups;
      m_grouping.groups(
        std::back_inserter(parallel_groups));

      Indices orthogonal_group;
      std::vector<bool> states(parallel_groups.size(), false);
      for (std::size_t i = 0; i < parallel_groups.size(); ++i) {
        if (states[i]) continue;

        CGAL_assertion(parallel_groups[i].size() > 0);
        const std::size_t si_index = parallel_groups[i][0];
        const auto& si = get(
          m_segment_map, *(m_input_range.begin() + si_index));

        states[i] = true;
        orthogonal_group.clear();
        for (const std::size_t seg_index : parallel_groups[i]) {
          orthogonal_group.push_back(seg_index);
        }

        traverse_group(
          preserve_order, i, si, parallel_groups,
          states, orthogonal_group);
        m_orthogonal_groups.push_back(orthogonal_group);
      }
      CGAL_assertion(
        m_orthogonal_groups.size() <= parallel_groups.size());
    }

    void traverse_group(
      const bool preserve_order,
      const std::size_t i,
      const Segment_2& si,
      const std::vector<Indices>& parallel_groups,
      std::vector<bool>& states,
      Indices& orthogonal_group) const {

      for (std::size_t j = i + 1; j < parallel_groups.size(); ++j) {
        if (states[j]) continue;

        CGAL_assertion(parallel_groups[j].size() > 0);
        const std::size_t sj_index = parallel_groups[j][0];
        const auto& sj = get(
          m_segment_map, *(m_input_range.begin() + sj_index));

        const FT angle_2 = internal::angle_2(si, sj);
        if (angle_2 <= m_max_angle ||
            angle_2 >= FT(90) - m_max_angle) {

          states[j] = true;
          for (const std::size_t seg_index : parallel_groups[j]) {
            orthogonal_group.push_back(seg_index);
          }
        } else if (preserve_order) return;
      }
    }
  };

} // namespace internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ORTHOGONAL_GROUPS_2_H
