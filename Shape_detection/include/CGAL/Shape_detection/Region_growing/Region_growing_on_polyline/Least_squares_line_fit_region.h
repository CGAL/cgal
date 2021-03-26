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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_LEAST_SQUARES_LINE_FIT_REGION_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_LEAST_SQUARES_LINE_FIT_REGION_H

#include <CGAL/license/Shape_detection.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/region_growing_traits.h>

namespace CGAL {
namespace Shape_detection {
namespace Polyline {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Least_squares_line_fit_region {

  public:
    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Point_type = typename Point_map::value_type;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// @}

  private:
    using Polyline_traits = typename std::conditional<
      std::is_same<typename Traits::Point_2, Point_type>::value,
      internal::Polyline_traits_2<Traits>,
      internal::Polyline_traits_3<Traits> >::type;

    using Point = typename Polyline_traits::Point;
    using Vector = typename Polyline_traits::Vector;
    using Line = typename Polyline_traits::Line;

    using Squared_length = typename Polyline_traits::Compute_squared_length;
    using Squared_distance = typename Polyline_traits::Compute_squared_distance;
    using Scalar_product = typename Polyline_traits::Compute_scalar_product;

  public:
    /// \name Initialization
    /// @{

    template<typename NamedParameters>
    Least_squares_line_fit_region(
      const InputRange& input_range,
      const NamedParameters& np,
      const PointMap point_map = PointMap(),
      const GeomTraits traits = GeomTraits()) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_polyline_traits(traits),
    m_squared_length(m_polyline_traits.compute_squared_length_object()),
    m_squared_distance(m_polyline_traits.compute_squared_distance_object()),
    m_scalar_product(m_polyline_traits.compute_scalar_product_object()) {

      CGAL_precondition(input_range.size() > 0);
      m_distance_threshold = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::distance_threshold), FT(1));
      CGAL_precondition(m_distance_threshold >= FT(0));

      const FT angle_deg_threshold = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::angle_deg_threshold), FT(25));
      CGAL_precondition(angle_deg_threshold >= FT(0) && angle_deg_threshold <= FT(90));

      m_min_region_size = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::min_region_size), 2);
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
      const std::size_t index1, const std::size_t index2,
      const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      CGAL_precondition(index1 < m_input_range.size());
      CGAL_precondition(index2 < m_input_range.size());
      const auto& key1 = *(m_input_range.begin() + index1);
      const auto& key2 = *(m_input_range.begin() + index2);
      const Point& input_point = get(m_point_map, key1);
      const Point& query_point = get(m_point_map, key2);

      if (region.size() == 1) { // update new reference line and direction
        CGAL_precondition(input_point != query_point);
        m_line_of_best_fit = Line(input_point, query_point);
        m_direction_of_best_fit = m_line_of_best_fit.to_vector();
        return true;
      }

      CGAL_precondition(region.size() >= 2);
      if (input_point == query_point) return true;
      CGAL_precondition(input_point != query_point);
      const Vector query_direction(input_point, query_point);

      CGAL_precondition(m_line_of_best_fit != Line());
      CGAL_precondition(m_direction_of_best_fit != Vector());

      const FT squared_distance_to_fitted_line =
        m_squared_distance(query_point, m_line_of_best_fit);
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
        m_line_of_best_fit = Line();
        m_direction_of_best_fit = Vector();
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

      // The best fit line will be a line fitted to all region points with
      // its direction being the line's direction.
      CGAL_precondition(region.size() > 0);
      const Line line_of_best_fit =
        m_polyline_traits.create_line_from_points(
          m_input_range, m_point_map, region).first;
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

    const Point_map m_point_map;
    const Polyline_traits m_polyline_traits;

    const Squared_length m_squared_length;
    const Squared_distance m_squared_distance;
    const Scalar_product m_scalar_product;

    Line m_line_of_best_fit;
    Vector m_direction_of_best_fit;
  };

} // namespace Polyline
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_LEAST_SQUARES_LINE_FIT_REGION_H
