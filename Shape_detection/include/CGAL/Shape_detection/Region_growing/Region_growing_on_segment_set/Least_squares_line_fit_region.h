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
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LEAST_SQUARES_LINE_FIT_REGION_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LEAST_SQUARES_LINE_FIT_REGION_H

#include <CGAL/license/Shape_detection.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/region_growing_traits.h>

namespace CGAL {
namespace Shape_detection {
namespace Segment_set {

  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap>
  class Least_squares_line_fit_region {

  public:
    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    using Segment_type = typename Segment_map::value_type;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// @}

  private:
    using Segment_set_traits = typename std::conditional<
      std::is_same<typename Traits::Segment_2, Segment_type>::value,
      internal::Region_growing_traits_2<Traits>,
      internal::Region_growing_traits_3<Traits> >::type;

    using Point = typename Segment_set_traits::Point;
    using Segment = typename Segment_set_traits::Segment;
    using Vector = typename Segment_set_traits::Vector;
    using Line = typename Segment_set_traits::Line;

    using Squared_length = typename Segment_set_traits::Compute_squared_length;
    using Squared_distance = typename Segment_set_traits::Compute_squared_distance;
    using Scalar_product = typename Segment_set_traits::Compute_scalar_product;

  public:
    /// \name Initialization
    /// @{

    template<typename NamedParameters>
    Least_squares_line_fit_region(
      const InputRange& input_range,
      const NamedParameters& np,
      const SegmentMap segment_map = SegmentMap(),
      const GeomTraits traits = GeomTraits()) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_segment_set_traits(traits),
    m_squared_length(m_segment_set_traits.compute_squared_length_object()),
    m_squared_distance(m_segment_set_traits.compute_squared_distance_object()),
    m_scalar_product(m_segment_set_traits.compute_scalar_product_object()) {

      CGAL_precondition(input_range.size() > 0);
      m_distance_threshold = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::distance_threshold), FT(1));
      CGAL_precondition(m_distance_threshold >= FT(0));

      const FT angle_deg_threshold = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::angle_deg_threshold), FT(25));
      CGAL_precondition(angle_deg_threshold >= FT(0) && angle_deg_threshold <= FT(90));

      m_min_region_size = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::min_region_size), 1);
      CGAL_precondition(m_min_region_size > 0);

      const FT cos_value_threshold = static_cast<FT>(std::cos(CGAL::to_double(
        (angle_deg_threshold * static_cast<FT>(CGAL_PI)) / FT(180))));
      m_cos_value_threshold = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::cos_value_threshold), cos_value_threshold);
      CGAL_precondition(m_cos_value_threshold >= FT(0) && m_cos_value_threshold <= FT(1));
    }

    /// @}

    /// \name Access
    /// @{

    bool is_part_of_region(
      const std::size_t,
      const std::size_t query_index,
      const std::vector<std::size_t>&) const {

      CGAL_precondition(query_index < m_input_range.size());
      const auto& key = *(m_input_range.begin() + query_index);
      const Segment& query_segment = get(m_segment_map, key);
      const Point& query_source = query_segment.source();
      const Point& query_target = query_segment.target();
      CGAL_precondition(query_source != query_target);
      const Vector query_direction(query_source, query_target);
      CGAL_precondition(query_direction != Vector());

      const FT squared_distance_to_fitted_line =
        get_max_squared_distance(query_segment);
      const FT squared_distance_threshold =
        m_distance_threshold * m_distance_threshold;

      const FT cos_value =
        m_scalar_product(query_direction, m_direction_of_best_fit);
      const FT squared_cos_value = cos_value * cos_value;

      FT squared_cos_value_threshold =
        m_cos_value_threshold * m_cos_value_threshold;
      squared_cos_value_threshold *= m_squared_length(query_direction);
      squared_cos_value_threshold *= m_squared_length(m_direction_of_best_fit);

      return (
        ( squared_distance_to_fitted_line <= squared_distance_threshold ) &&
        ( squared_cos_value >= squared_cos_value_threshold ));
    }

    inline bool is_valid_region(const std::vector<std::size_t>& region) const {
      return (region.size() >= m_min_region_size);
    }

    void update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference line and direction
        const std::size_t segment_index = region[0];
        CGAL_precondition(segment_index < m_input_range.size());

        // The best fit line will be a line obtained from this segment
        // with the same direction.
        const auto& key = *(m_input_range.begin() + segment_index);
        const Segment& segment = get(m_segment_map, key);
        const Point& source = segment.source();
        const Point& target = segment.target();
        CGAL_precondition(source != target);

        m_line_of_best_fit = Line(source, target);
        m_direction_of_best_fit = m_line_of_best_fit.to_vector();

      } else { // update reference line and direction
        CGAL_precondition(region.size() >= 2);
        std::tie(m_line_of_best_fit, m_direction_of_best_fit) =
          get_line_and_direction(region);
      }
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    std::pair<Line, Vector> get_line_and_direction(
      const std::vector<std::size_t>& region) const {

      // The best fit line will be a line fitted to all region segments with
      // its direction being the line's direction.
      CGAL_precondition(region.size() > 0);
      const Line line_of_best_fit =
        m_segment_set_traits.create_line(
          m_input_range, m_segment_map, region).first;
      const Vector direction_of_best_fit =
        line_of_best_fit.to_vector();

      return std::make_pair(line_of_best_fit, direction_of_best_fit);
    }
    /// \endcond

  private:
    const Input_range& m_input_range;

    FT m_distance_threshold;
    FT m_cos_value_threshold;
    std::size_t m_min_region_size;

    const Segment_map m_segment_map;
    const Segment_set_traits m_segment_set_traits;

    const Squared_length m_squared_length;
    const Squared_distance m_squared_distance;
    const Scalar_product m_scalar_product;

    Line m_line_of_best_fit;
    Vector m_direction_of_best_fit;

    // The maximum squared distance from the vertices of the segment
    // to the best fit line.
    FT get_max_squared_distance(const Segment& segment) const {

      const Point& source = segment.source();
      const Point& target = segment.target();
      const FT squared_distance_source =
        m_squared_distance(source, m_line_of_best_fit);
      const FT squared_distance_target =
        m_squared_distance(target, m_line_of_best_fit);
      return (CGAL::max)(squared_distance_source, squared_distance_target);
    }
  };

} // namespace Segment_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LEAST_SQUARES_LINE_FIT_REGION_H
