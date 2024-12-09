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

#ifndef CGAL_SHAPE_REGULARIZATION_UNIQUE_SEGMENTS_2_H
#define CGAL_SHAPE_REGULARIZATION_UNIQUE_SEGMENTS_2_H

#include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/internal/Collinear_groups_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap>
  class Unique_segments_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;

    using FT = typename Traits::FT;
    using Line_2 = typename Traits::Line_2;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Segment_2 = typename Traits::Segment_2;
    using Indices = std::vector<std::size_t>;

    using CGroups_2 = Collinear_groups_2<Traits, Input_range, Segment_map>;

    template<typename NamedParameters>
    Unique_segments_2(
      const InputRange& input_range,
      const NamedParameters& np,
      const SegmentMap segment_map,
      const GeomTraits&) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_grouping(
      input_range, np, segment_map, GeomTraits()) {

      make_unique_segments();
    }

    template<typename OutputIterator>
    OutputIterator segments(OutputIterator segments) const {
      for (const auto& segment : m_segments) {
        *(segments++) = segment;
      }
      return segments;
    }

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;
    const CGroups_2 m_grouping;

    std::vector<Segment_2> m_segments;

    void make_unique_segments() {

      std::vector<Indices> collinear_groups;
      m_grouping.groups(
        std::back_inserter(collinear_groups));
      m_segments.reserve(collinear_groups.size());

      for (const auto& collinear_group : collinear_groups) {
        handle_collinear_group(collinear_group);
      }
      CGAL_assertion(m_segments.size() == collinear_groups.size());
    }

    void handle_collinear_group(
      const Indices& collinear_group) {

      CGAL_assertion(collinear_group.size() > 0);
      if (collinear_group.size() == 1) {
        const std::size_t seg_index = collinear_group[0];
        CGAL_assertion(seg_index < m_input_range.size());
        const auto& first = get(m_segment_map,
          *(m_input_range.begin() + seg_index));
        m_segments.push_back(first);
        return;
      }
      m_segments.push_back(
        find_weighted_segment(collinear_group));
    }

    Segment_2 find_weighted_segment(
      const std::vector<std::size_t>& collinear_group) const {

      std::vector<FT> weights;
        compute_distance_weights(collinear_group, weights);

      const std::size_t longest =
        find_longest_segment(collinear_group);
      CGAL_assertion(longest < m_input_range.size());
      const auto& ref_segment = get(m_segment_map,
        *(m_input_range.begin() + longest));

      Segment_2 weighted = compute_weighted_segment(
        collinear_group, weights, ref_segment);
      if (weighted.source() == weighted.target()) {
        weighted = ref_segment;
      }

      const Vector_2 ref_vector = weighted.to_vector();
      const Line_2 ref_line = Line_2(weighted.source(), weighted.target());

      Point_2 source, target;
      compute_boundary_points(
        collinear_group, ref_vector, source, target);
      if (source != target) {

        source = ref_line.projection(source);
        target = ref_line.projection(target);
        weighted = Segment_2(source, target);
      }
      return weighted;
    }

    void compute_distance_weights(
      const std::vector<std::size_t>& collinear_group,
      std::vector<FT>& weights) const {

      CGAL_assertion(
        collinear_group.size() > 0);
      weights.clear();
      weights.reserve(collinear_group.size());

      FT sum_distance = FT(0);
      for (const std::size_t seg_index : collinear_group) {
        CGAL_assertion(seg_index < m_input_range.size());
        const auto& segment = get(m_segment_map,
          *(m_input_range.begin() + seg_index));

        const FT sq_distance = segment.squared_length();
        sum_distance += sq_distance;
        weights.push_back(sq_distance);
      }

      CGAL_assertion(sum_distance > FT(0));
      for (auto& weight : weights) {
        weight /= sum_distance;
      }
      CGAL_assertion(
        weights.size() == collinear_group.size());
    }

    std::size_t find_longest_segment(
      const Indices& collinear_group) const {

      FT max_length = -FT(1);
      std::size_t longest = std::size_t(-1);

      for (const std::size_t seg_index : collinear_group) {
        CGAL_assertion(seg_index < m_input_range.size());
        const auto& segment = get(m_segment_map,
          *(m_input_range.begin() + seg_index));

        const FT length = segment.squared_length();
        if (length > max_length) {
          longest = seg_index; max_length = length;
        }
      }
      CGAL_assertion(longest != std::size_t(-1));
      return longest;
    }

    Segment_2 compute_weighted_segment(
      const Indices& collinear_group,
      const std::vector<FT>& weights,
      const Segment_2& ref_segment) const {

      const auto& sref = ref_segment.source();
      const auto& tref = ref_segment.target();

      const auto center = CGAL::midpoint(sref, tref);
      CGAL_assertion(
        weights.size() == collinear_group.size());
      Vector_2 dir = Vector_2(FT(0), FT(0));
      for (std::size_t i = 0; i < weights.size(); ++i) {
        const FT weight = weights[i];

        const std::size_t seg_index = collinear_group[i];
        CGAL_assertion(seg_index < m_input_range.size());
        const auto& segment = get(m_segment_map,
          *(m_input_range.begin() + seg_index));

        const Line_2 line = Line_2(
          segment.source(), segment.target());
        const Point_2 proj = line.projection(center);

        const Vector_2 v = Vector_2(center, proj);
        dir += v * weight;
      }

      const Point_2 source = sref + dir;
      const Point_2 target = tref + dir;

      return Segment_2(source, target);
    }

    void compute_boundary_points(
      const std::vector<std::size_t>& collinear_group,
      const Vector_2& ref_vector,
      Point_2& p, Point_2& q) const {

      const FT max_value = internal::max_value<FT>();
      FT min_proj_value =  max_value;
      FT max_proj_value = -max_value;

      CGAL_assertion(collinear_group.size() > 0);
      CGAL_assertion(collinear_group[0] < m_input_range.size());
      const auto& first = get(m_segment_map,
        *(m_input_range.begin() + collinear_group[0]));

      const Point_2& ref_point = first.source();
      for (const std::size_t seg_index : collinear_group) {
        CGAL_assertion(seg_index < m_input_range.size());
        const auto& segment = get(m_segment_map,
          *(m_input_range.begin() + seg_index));

        const auto& source = segment.source();
        const auto& target = segment.target();

        recompute_end_points(
          source, ref_point, ref_vector,
          min_proj_value, max_proj_value, p, q);
        recompute_end_points(
          target, ref_point, ref_vector,
          min_proj_value, max_proj_value, p, q);
      }
    }

    void recompute_end_points(
      const Point_2& query,
      const Point_2& ref_point,
      const Vector_2& ref_vector,
      FT& min_proj_value,
      FT& max_proj_value,
      Point_2& p,
      Point_2& q) const {

      const Vector_2 curr_vector(ref_point, query);
      const FT value = CGAL::scalar_product(curr_vector, ref_vector);
      if (value < min_proj_value) {
        min_proj_value = value;
        p = query;
      }
      if (value > max_proj_value) {
        max_proj_value = value;
        q = query;
      }
    }
  };

} // namespace internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_UNIQUE_SEGMENTS_2_H
