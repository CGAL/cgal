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

#ifndef CGAL_SHAPE_REGULARIZATION_CLOSED_CONTOUR_2_H
#define CGAL_SHAPE_REGULARIZATION_CLOSED_CONTOUR_2_H

#include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Contour_base_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<
  typename ContourDirections,
  typename GeomTraits>
  class Closed_contour_2 {

  public:
    using Contour_directions = ContourDirections;
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;
    using Base = internal::Contour_base_2<Traits>;
    using Segment_wrapper_2 = typename Base::Segment_wrapper_2;

    Closed_contour_2(
      const Contour_directions& estimator,
      const FT max_offset_2) :
    m_estimator(estimator),
    m_max_offset_2(max_offset_2)
    { }

    template<
    typename Input_range,
    typename Point_map>
    void initialize(
      const Input_range& input_range, const Point_map point_map) {

      m_base.initialize_closed(
        input_range, point_map, m_wraps);
      CGAL_assertion(m_wraps.size() == input_range.size());
    }

    template<typename OutputIterator>
    OutputIterator regularize(OutputIterator contour) {

      if (m_wraps.size() < 3) return contour;
      CGAL_assertion(m_wraps.size() >= 3);

      rotate_contour(m_wraps);
      if (verbose()) {
        m_base.export_polylines(m_wraps, "rotated");
      }

      bool success = optimize_contour(m_wraps);
      if (!success) return contour;
      if (verbose()) {
        m_base.export_polylines(m_wraps, "optimized");
      }

      success = connect_contour(m_wraps);
      if (!success) return contour;
      if (verbose()) {
        m_base.export_polylines(m_wraps, "connected");
      }

      return update_input(m_wraps, contour);
    }

  private:
    const Contour_directions& m_estimator;
    const FT m_max_offset_2;
    const Base m_base;

    std::vector<Segment_wrapper_2> m_wraps;

    bool verbose() const {
      return m_base.verbose();
    }

    void rotate_contour(
      std::vector<Segment_wrapper_2>& wraps) const {

      for (std::size_t i = 0; i < wraps.size(); ++i) {
        auto& wrap = wraps[i];
        m_estimator.orient(i, wrap.segment);
      }
    }

    bool optimize_contour(
      std::vector<Segment_wrapper_2>& wraps) const {

      // Clean.
      m_base.remove_zero_length_segments(wraps);
      if (wraps.size() < 4) return false; // should be at least a quad
      CGAL_assertion(wraps.size() >= 4);

      // Merge parallel/collinear segments.
      std::vector<Segment_2> segments;
      m_base.create_unique_segments(m_max_offset_2, wraps, segments);
      if (verbose()) {
        std::cout <<
          "* number of segments (merging) = " << segments.size() << std::endl;
      }
      if (wraps.size() < 4) return false;

      // Add orthogonal segments.
      create_orthogonal_segments(segments, wraps);
      if (verbose()) {
        std::cout <<
          "* number of segments (orthogonal) = " << wraps.size() << std::endl;
      }
      if (wraps.size() < 4) return false;

      return true;
    }

    void create_orthogonal_segments(
      const std::vector<Segment_2>& segments,
      std::vector<Segment_wrapper_2>& wraps) const {

      Segment_wrapper_2 wrap;
      const std::size_t n = segments.size();

      wraps.clear();
      wraps.reserve(n);

      std::size_t count = 0;
      for (std::size_t i = 0; i < n; ++i) {
        const auto& segmenti = segments[i];

        wrap.index = count; ++count;
        wrap.segment = segmenti;
        wraps.push_back(wrap);

        const std::size_t j = (i + 1) % n;
        const auto& segmentj = segments[j];

        if (m_base.is_parallel_segment(segmenti, segmentj)) {
          wrap.index = count; ++count;
          m_base.create_average_orth(
            segmenti, segmentj, wrap.segment);
          wraps.push_back(wrap);
        }
      }
      CGAL_assertion(wraps.size() >= segments.size());
    }

    bool connect_contour(
      std::vector<Segment_wrapper_2>& wraps) const {

      CGAL_assertion(wraps.size() >= 4);
      if (wraps.size() < 4) return false;

      intersect_segments(wraps);
      if (wraps.size() < 4) return false;

      // Experimental code.
      // make_segments_collinear(wraps);
      // intersect_segments(wraps);
      // if (wraps.size() < 4) return false;

      return true;
    }

    void intersect_segments(
      std::vector<Segment_wrapper_2>& wraps) const {

      const std::size_t n = wraps.size();
      for (std::size_t i = 0; i < n; ++i) {

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        auto& si = wraps[i].segment;
        const auto& sm = wraps[im].segment;
        const auto& sp = wraps[ip].segment;

        m_base.intersect_segment(sm, si, sp);
      }
    }

    // Experimental function.
    void make_segments_collinear(
      std::vector<Segment_wrapper_2>& wraps) const {

      Segment_wrapper_2 wrap;
      std::vector<Segment_wrapper_2> clean, group;

      const std::size_t n = wraps.size();
      const FT sq_max_length = m_max_offset_2 * m_max_offset_2;
      for (std::size_t i = 0; i < n; ++i) {
        group.push_back(wraps[i]);

        const std::size_t j = (i + 1) % n;
        const FT sq_length = wraps[j].segment.squared_length();
        if (sq_length < sq_max_length) {
          ++i; continue;
        } else {
          m_base.parallel_segments_to_segment(group, wrap.segment);
          clean.push_back(wrap); group.clear();
        }
      }

      const std::size_t before = wraps.size();
      wraps = clean;
      const std::size_t after = wraps.size();
      CGAL_assertion(after <= before);

      if (verbose()) {
        std::cout <<
          "* segments before/after = " << before << "/" << after << std::endl;
      }
    }

    template<typename OutputIterator>
    OutputIterator update_input(
      const std::vector<Segment_wrapper_2>& wraps,
      OutputIterator contour) const {

      for (std::size_t i = 0; i < wraps.size(); ++i) {
        const auto& wrap = wraps[i];
        *(++contour) = wrap.segment.source();
      }
      return contour;
    }
  };

} // internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CLOSED_CONTOUR_2_H
