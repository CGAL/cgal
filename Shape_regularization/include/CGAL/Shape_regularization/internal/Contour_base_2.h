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

#ifndef CGAL_SHAPE_REGULARIZATION_CONTOUR_BASE_2_H
#define CGAL_SHAPE_REGULARIZATION_CONTOUR_BASE_2_H

#include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/internal/Segment_wrapper_2.h>
#include <CGAL/Shape_regularization/internal/Unique_segments_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<typename GeomTraits>
  class Contour_base_2 {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_2 = typename Traits::Vector_2;
    using Segment_2 = typename Traits::Segment_2;
    using Direction_2 = typename Traits::Direction_2;
    using Line_2 = typename Traits::Line_2;
    using Intersect_2 = typename Traits::Intersect_2;

    using FT_pair = std::pair<FT, FT>;
    using Segment_wrapper_2 = internal::Segment_wrapper_2<Traits>;
    using Segment_wrappers_2 = std::vector<Segment_wrapper_2>;
    using Polyline = std::vector<Point_3>;

    Contour_base_2() :
    m_verbose(false),
    m_angle_threshold_2(FT(5))
    { }

    /////////////////////
    // Debug and Save. //
    /////////////////////

    void export_polylines(
      const std::vector<Segment_wrapper_2>& wraps,
      const std::string path) const {

      std::vector<Segment_2> segments;
      segments.reserve(wraps.size());
      for (const auto& wrap : wraps)
        segments.push_back(wrap.segment);
      export_polylines(segments, path);
    }

    void export_polylines(
      const std::vector<Segment_2>& segments,
      const std::string path) const {

      std::vector<Polyline> polylines(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const auto& s = segments[i].source();
        const auto& t = segments[i].target();

        polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
        polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
      }
      export_polylines(polylines, path);
    }

    void export_polylines(
      const std::vector<Polyline>& polylines,
      const std::string path) const {

      if (polylines.size() == 0) return;
      std::stringstream out;
      out.precision(20);

      for (std::size_t i = 0; i < polylines.size(); ++i) {
        const auto& polyline = polylines[i];

        out << polyline.size() << " ";
        for (std::size_t j = 0; j < polyline.size(); ++j) {
          out << polyline[j] << " ";
        }
        out << std::endl;
      }
      save(out, path + ".polylines");
    }

    void save(
      const std::stringstream& out,
      const std::string path) const {

      std::ofstream file(path.c_str(), std::ios_base::out);
      CGAL::IO::set_ascii_mode(file);
      if (!file) {
        std::cout <<
          "Error: cannot save the file: " << path << std::endl;
        return;
      }
      file << out.str() << std::endl; file.close();
      std::cout <<
        "* segments are saved in " << path << std::endl;
    }

    ////////////////////////
    // General Utilities. //
    ////////////////////////

    bool verbose() const {
      return m_verbose;
    }

    const FT get_angle_threshold_2() const {
      return m_angle_threshold_2;
    }

    void sort_segments_by_length(
      const std::vector<Segment_wrapper_2>& wraps,
      std::vector<std::size_t>& sorted) const {

      sorted.clear();
      sorted.reserve(wraps.size());
      for (std::size_t i = 0; i < wraps.size(); ++i) {
        sorted.push_back(i);
      }

      std::sort(sorted.begin(), sorted.end(),
      [&wraps](const std::size_t i, const std::size_t j) -> bool {

        const FT length_1 = wraps[i].segment.squared_length();
        const FT length_2 = wraps[j].segment.squared_length();
        return length_1 > length_2;
      });
    }

    /////////////////////
    // Initialization. //
    /////////////////////

    template<
    typename Input_range,
    typename Point_map>
    void initialize_closed(
      const Input_range& input_range,
      const Point_map point_map,
      std::vector<Segment_wrapper_2>& wraps) const {

      CGAL_precondition(input_range.size() >= 3);
      const std::size_t n = input_range.size();

      wraps.clear();
      wraps.reserve(n);

      Segment_wrapper_2 wrap;
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ip = (i + 1) % n;

        const auto& source = get(point_map, *(input_range.begin() + i));
        const auto& target = get(point_map, *(input_range.begin() + ip));

        wrap.index = i;
        wrap.segment = Segment_2(source, target);
        wrap.set_direction(wrap.segment);
        wraps.push_back(wrap);
      }
      CGAL_assertion(wraps.size() == n);
    }

    template<
    typename Input_range,
    typename Point_map>
    void initialize_open(
      const Input_range& input_range,
      const Point_map point_map,
      std::vector<Segment_wrapper_2>& wraps) const {

      CGAL_precondition(input_range.size() >= 2);
      const std::size_t n = input_range.size();

      wraps.clear();
      wraps.reserve(n);

      Segment_wrapper_2 wrap;
      for (std::size_t i = 0; i < n - 1; ++i) {
        const std::size_t ip = i + 1;

        const auto& source = get(point_map, *(input_range.begin() + i));
        const auto& target = get(point_map, *(input_range.begin() + ip));

        wrap.index = i;
        wrap.segment = Segment_2(source, target);
        wrap.set_direction(wrap.segment);
        wraps.push_back(wrap);
      }
      CGAL_assertion(wraps.size() == n - 1);
    }

    /////////////////
    // Directions. //
    /////////////////

    void unify_along_contours_closed(
      std::vector<Segment_wrapper_2>& wraps,
      std::vector<std::size_t>& assigned) const {

      const std::size_t n = wraps.size();
      CGAL_assertion(assigned.size() == n);
      for (std::size_t i = 0; i < n; ++i) {
        auto& wrap = wraps[i];
        if (wrap.is_used) continue;

        std::size_t im = (i + n - 1) % n;
        std::size_t ip = (i + 1) % n;

        bool stop = false;
        std::size_t max_count = 0;
        do {

          if (wraps[im].is_used) {
            assigned[i] = assigned[im];
            wrap.is_used = true;
            break;
          }

          if (wraps[ip].is_used) {
            assigned[i] = assigned[ip];
            wrap.is_used = true;
            break;
          }

          im = (im + n - 1) % n;
          ip = (ip + 1) % n;

          if (im == i || ip == i) {
            stop = true;
          }
          ++max_count;

        } while (!stop && max_count < n);
        if (stop || max_count >= n) {
          std::cerr <<
            "Warning: revert back to the first direction!" << std::endl;
          assigned[i] = 0;
        }
      }
    }

    void correct_directions_closed(
      const std::vector<Segment_wrapper_2>& wraps,
      std::vector<std::size_t>& assigned) const {

      const std::size_t n = wraps.size();
      std::vector<std::size_t> clean;
      clean.reserve(n);

      CGAL_assertion(assigned.size() == n);
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        const std::size_t dm = assigned[im];
        const std::size_t di = assigned[i];
        const std::size_t dp = assigned[ip];

        if (dm != std::size_t(-1) && dm == dp && di != dm) {
          clean.push_back(dm);
        } else {
          clean.push_back(di);
        }
      }
      assigned = clean;
    }

    void unify_along_contours_open(
      std::vector<Segment_wrapper_2>& wraps,
      std::vector<std::size_t>& assigned) const {

      CGAL_assertion(assigned.size() == wraps.size());
      const int n = static_cast<int>(wraps.size());

      for (int i = 0; i < n; ++i) {
        auto& wrap = wraps[i];
        if (wrap.is_used) continue;

        int im = -1;
        if (i > 0) im = i - 1;
        int ip = -1;
        if (i < n - 1) ip = i + 1;

        bool stop = false;
        int max_count = 0;
        do {

          if (im != -1 && wraps[im].is_used) {
            CGAL_assertion(i >= 0 && i < n);
            CGAL_assertion(im >= 0 && im < n);
            assigned[i] = assigned[im];
            wrap.is_used = true;
            break;
          }

          if (ip != -1 && wraps[ip].is_used) {
            CGAL_assertion(i >= 0 && i < n);
            CGAL_assertion(ip >= 0 && ip < n);
            assigned[i] = assigned[ip];
            wrap.is_used = true;
            break;
          }

          if (stop) break;
          if (im != -1 && im > 0) {
            im = im - 1;
          }
          if (ip != -1 && ip < n - 1) {
            ip = ip + 1;
          }

          if (im == 0 || ip == n - 1) {
            stop = true;
          }
          ++max_count;

        } while (max_count < n);
        if (stop || max_count >= n) {
          std::cerr <<
            "Warning: revert back to the first direction!" << std::endl;
          assigned[i] = 0;
        }
      }
    }

    void correct_directions_open(
      std::vector<Segment_wrapper_2>& wraps,
      std::vector<std::size_t>& assigned) const {

      CGAL_assertion(assigned.size() == wraps.size());
      const int n = static_cast<int>(wraps.size());
      std::vector<std::size_t> clean;
      clean.reserve(n);

      for (int i = 0; i < n; ++i) {

        if (i == 0) {
          const int ip = 1;
          CGAL_assertion(ip >= 0 && ip < n);
          const std::size_t di = assigned[i];
          const std::size_t dp = assigned[ip];
          if (di != dp) {
            clean.push_back(dp);
          } else {
            clean.push_back(di);
          }
          continue;
        }

        if (i == n - 1) {
          const int im = n - 2;
          CGAL_assertion(im >= 0 && im < n);
          const std::size_t dm = assigned[im];
          const std::size_t di = assigned[i];
          if (di != dm) {
            clean.push_back(dm);
          } else {
            clean.push_back(di);
          }
          continue;
        }

        const int im = i - 1;
        const int ip = i + 1;
        CGAL_assertion(im >= 0 && im < n);
        const std::size_t dm = assigned[im];
        CGAL_assertion(i >= 0 && i < n);
        const std::size_t di = assigned[i];
        CGAL_assertion(ip >= 0 && ip < n);
        const std::size_t dp = assigned[ip];

        if (dm != std::size_t(-1) && dm == dp && di != dm) {
          clean.push_back(dm);
        } else {
          clean.push_back(di);
        }
      }
      assigned = clean;
    }

    void readjust_directions(
      const std::vector<Segment_wrapper_2>& wraps,
      const std::vector<std::size_t>& assigned,
      std::vector<Direction_2>& directions) const {

      std::vector<FT> angles, counts;
      create_average_angles(wraps, assigned, directions,
      angles, counts);

      CGAL_assertion(angles.size() == counts.size());
      CGAL_assertion(angles.size() == directions.size());

      for (std::size_t k = 0; k < angles.size(); ++k) {
        if (counts[k] == FT(0)) continue;
        angles[k] /= counts[k];

        const FT angle_deg = angles[k];
        internal::rotate_direction_2(angle_deg, directions[k]);
      }
    }

    void create_average_angles(
      const std::vector<Segment_wrapper_2>& wraps,
      const std::vector<std::size_t>& assigned,
      const std::vector<Direction_2>& directions,
      std::vector<FT>& angles,
      std::vector<FT>& counts) const {

      CGAL_assertion(directions.size() > 0);

      angles.clear();
      angles.resize(directions.size(), FT(0));

      counts.clear();
      counts.resize(directions.size(), FT(0));

      for (std::size_t i = 0; i < wraps.size(); ++i) {
        const auto& wrap = wraps[i];
        if (!wrap.is_valid_direction) continue;

        const std::size_t direction_index = assigned[i];
        CGAL_assertion(direction_index != std::size_t(-1));

        // Experimental code.
        // const auto& di = directions[direction_index];
        // const auto& dj = wrap.direction;

        const auto& di = wrap.direction;
        const auto& dj = directions[direction_index];
        const FT angle = internal::mod90_angle_2(di, dj);

        angles[direction_index] += angle;
        counts[direction_index] += FT(1);
      }
    }

    void apply_rotation_to_segment(
      const std::vector<FT_pair>& bounds,
      const std::vector<Direction_2>& directions,
      const std::vector<std::size_t>& assigned,
      const std::size_t query_index,
      Segment_2& segment) const {

      CGAL_assertion(assigned.size() > 0);
      CGAL_assertion(bounds.size() == directions.size());
      CGAL_assertion(query_index < assigned.size());

      const std::size_t direction_index = assigned[query_index];
      if (direction_index == std::size_t(-1)) {
        return;
      }

      CGAL_assertion(direction_index < directions.size());
      const auto& ref_direction = directions[direction_index];
      const auto& ref_bounds = bounds[direction_index];

      auto v = segment.to_vector();
      const Direction_2 seg_direction =
        internal::direction_2(v);
      rotate_segment(
        ref_bounds, ref_direction, seg_direction, segment);
    }

    void rotate_segment(
      const FT_pair& /* bounds */,
      const Direction_2& ref_direction,
      const Direction_2& seg_direction,
      Segment_2& segment) const {

      // Experimental code.
      // const FT angle_deg = internal::compute_angle_2(
      //   ref_direction, seg_direction);
      // const FT converted = CGAL::abs(convert_angle_2(angle_deg));
      // if (converted <= bounds.first)
      //   internal::rotate_segment_2(
      //     angle_deg, FT(180), segment); // parallel case
      // if (converted >= bounds.second)
      //   internal::rotate_segment_2(
      //     angle_deg, FT(90), segment); // orthogonal case

      const FT angle_deg = internal::mod90_angle_2(
        seg_direction, ref_direction);
      internal::rotate_segment_2(
        angle_deg, FT(0), segment);
    }

    ///////////////
    // Contours. //
    ///////////////

    void remove_zero_length_segments(
      std::vector<Segment_wrapper_2>& wraps) const {

      std::vector<Segment_wrapper_2> clean;
      for (const auto& wrap : wraps) {
        if (wrap.segment.squared_length() > internal::tolerance<FT>()) {
          clean.push_back(wrap);
        }
      }
      wraps = clean;
    }

    void create_unique_segments(
      const FT max_offset_2,
      const std::vector<Segment_wrapper_2>& wraps,
      std::vector<Segment_2>& segments) const {

      using SegmentRange = std::vector<Segment_wrapper_2>;
      using SegmentMap = Wrap_segment_map<Traits>;
      using Unique_segments_2 = internal::Unique_segments_2<
        Traits, SegmentRange, SegmentMap>;

      const SegmentMap segment_map;
      const Unique_segments_2 unique(
        wraps, CGAL::parameters::
        maximum_angle(get_angle_threshold_2()).
        maximum_offset(max_offset_2).
        preserve_order(true),
        segment_map, Traits());

      segments.clear();
      unique.segments(
        std::back_inserter(segments));
    }

    std::pair<bool, bool> is_parallel_segment(
      const Segment_2& sm,
      const Segment_2& si,
      const Segment_2& sp) const {

      const FT angle_mi_2 = internal::angle_2(sm, si);
      const FT angle_pi_2 = internal::angle_2(si, sp);

      const bool source_cond = ( angle_mi_2 <= m_angle_threshold_2 );
      const bool target_cond = ( angle_pi_2 <= m_angle_threshold_2 );

      return std::make_pair(source_cond, target_cond);
    }

    bool is_parallel_segment(
      const Segment_2& si, const Segment_2& sp) const {

      const FT angle_pi_2 = internal::angle_2(si, sp);
      const bool target_cond = ( angle_pi_2 <= m_angle_threshold_2 );
      return target_cond;
    }

    void parallel_segments_to_segment(
      const std::vector<Segment_wrapper_2>& wraps,
      Segment_2& result) const {

      Segment_2 ref_segment = find_weighted_segment(wraps);
      const Line_2 line = Line_2(
        ref_segment.source(), ref_segment.target());

      std::vector<Point_2> points;
      for (const auto& wrap : wraps) {

        const Point_2 source = line.projection(wrap.segment.source());
        const Point_2 target = line.projection(wrap.segment.target());

        points.push_back(source);
        points.push_back(target);
      }
      update_segment(points, ref_segment);
      result = ref_segment;
    }

    Segment_2 find_weighted_segment(
      const std::vector<Segment_wrapper_2>& wraps) const {

      std::vector<FT> weights;
        compute_distance_weights(wraps, weights);
      const Segment_2 ref_segment =
        find_longest_segment(wraps);
      const Segment_2 weighted =
        compute_weighted_segment(wraps, weights, ref_segment);
      if (weighted.source() == weighted.target()) {
        return ref_segment;
      }
      return weighted;
    }

    void compute_distance_weights(
      const std::vector<Segment_wrapper_2>& wraps,
      std::vector<FT>& weights) const {

      CGAL_assertion(wraps.size() > 0);
      weights.clear();
      weights.reserve(wraps.size());

      FT sum_distance = FT(0);
      for (const auto& wrap : wraps) {
        const FT sq_distance = wrap.segment.squared_length();
        sum_distance += sq_distance;
        weights.push_back(sq_distance);
      }

      CGAL_assertion(sum_distance > FT(0));
      for (auto& weight : weights) {
        weight /= sum_distance;
      }
      CGAL_assertion(weights.size() == wraps.size());
    }

    Segment_2 find_longest_segment(
      const std::vector<Segment_wrapper_2>& wraps) const {

      FT max_length = -FT(1);
      std::size_t longest = std::size_t(-1);

      for (std::size_t i = 0; i < wraps.size(); ++i) {
        const auto& wrap = wraps[i];
        const FT length = wrap.segment.squared_length();
        if (length > max_length) {
          longest = i; max_length = length;
        }
      }
      CGAL_assertion(longest != std::size_t(-1));
      return wraps[longest].segment;
    }

    Segment_2 compute_weighted_segment(
      const std::vector<Segment_wrapper_2>& wraps,
      const std::vector<FT>& weights,
      const Segment_2& ref_segment) const {

      const auto& sref = ref_segment.source();
      const auto& tref = ref_segment.target();

      const auto center = CGAL::midpoint(sref, tref);
      CGAL_assertion(weights.size() == wraps.size());
      Vector_2 dir = Vector_2(FT(0), FT(0));
      for (std::size_t i = 0; i < weights.size(); ++i) {
        const FT weight = weights[i];

        const auto& wrap = wraps[i];
        const Line_2 line = Line_2(
          wrap.segment.source(), wrap.segment.target());
        const Point_2 proj = line.projection(center);

        const Vector_2 v = Vector_2(center, proj);
        dir += v * weight;
      }

      const Point_2 source = sref + dir;
      const Point_2 target = tref + dir;

      return Segment_2(source, target);
    }

    void update_segment(
      const std::vector<Point_2>& points,
      Segment_2& segment) const {

      FT min_proj_value = +internal::max_value<FT>();
      FT max_proj_value = -internal::max_value<FT>();

      const Vector_2 ref_vector = segment.to_vector();
      const Point_2 ref_point = internal::barycenter_2(points);

      Point_2 source, target;
      for (const auto& point : points) {
        const Vector_2 curr_vector(ref_point, point);
        const FT value = CGAL::scalar_product(curr_vector, ref_vector);

        if (value < min_proj_value) {
          min_proj_value = value;
          source = point;
        }
        if (value > max_proj_value) {
          max_proj_value = value;
          target = point;
        }
      }
      segment = Segment_2(source, target);
    }

    void intersect_segment(
      const Segment_2& sm,
      Segment_2& si) const {

      Point_2 source = si.source();
      Point_2 target = si.target();

      const Line_2 line_1 = Line_2(sm.source(), sm.target());
      const Line_2 line_2 = Line_2(si.source(), si.target());
      const bool success = intersect_2(line_1, line_2, source);

      if (!success) source = si.source();
      si = Segment_2(source, target);
    }

    void intersect_segment(
      const Segment_2& sm,
      Segment_2& si,
      const Segment_2& sp) const {

      Point_2 source = si.source();
      Point_2 target = si.target();

      const Line_2 line_1 = Line_2(sm.source(), sm.target());
      const Line_2 line_2 = Line_2(si.source(), si.target());
      const Line_2 line_3 = Line_2(sp.source(), sp.target());

      const bool success1 = intersect_2(line_1, line_2, source);
      const bool success2 = intersect_2(line_2, line_3, target);

      if (!success1) source = si.source();
      if (!success2) target = si.target();

      si = Segment_2(source, target);
    }

    void intersect_segment(
      Segment_2& si,
      const Segment_2& sp) const {

      Point_2 source = si.source();
      Point_2 target = si.target();

      const Line_2 line_1 = Line_2(si.source(), si.target());
      const Line_2 line_2 = Line_2(sp.source(), sp.target());
      const bool success = intersect_2(line_1, line_2, target);

      if (!success) target = si.target();
      si = Segment_2(source, target);
    }

    bool intersect_2(
      const Line_2& line_1,
      const Line_2& line_2,
      Point_2& in_point) const {

      typename CGAL::cpp11::result_of<Intersect_2(Line_2, Line_2)>::type result
      = CGAL::intersection(line_1, line_2);
      if (result) {
        if (const Line_2* line = std::get_if<Line_2>(&*result)) {
          return false;
        } else {
          const Point_2* point = std::get_if<Point_2>(&*result);
          in_point = *point; return true;
        }
      }
      return false;
    }

    void create_average_orth(
      const Segment_2& segmenti,
      const Segment_2& segmentj,
      Segment_2& orth) const {

      const Line_2 linei = Line_2(
        segmenti.source(), segmenti.target());
      const auto p = linei.projection(segmentj.source());
      const auto source = CGAL::midpoint(p, segmenti.target());

      const Line_2 linej = Line_2(
        segmentj.source(), segmentj.target());
      const auto q = linej.projection(segmenti.target());
      const auto target = CGAL::midpoint(q, segmentj.source());

      orth = Segment_2(source, target);
    }

  private:
    const bool m_verbose;
    const FT m_angle_threshold_2;
  };

} // internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CONTOUR_BASE_2_H
