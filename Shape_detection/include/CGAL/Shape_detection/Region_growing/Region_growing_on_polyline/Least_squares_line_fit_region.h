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

  /*!
    \ingroup PkgShapeDetectionRGOnPolyline

    \brief Region type based on the quality of the least squares line
    fit applied to polyline vertices.

    This class fits a line, using \ref PkgPrincipalComponentAnalysisDRef "PCA",
    to chunks of polyline vertices and controls the quality of this fit.
    If all quality conditions are satisfied, the chunk is accepted as a valid region,
    otherwise rejected.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam PointMap
    a model of `LValuePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Point_2` or `Kernel::Point_3`

    \cgalModels `RegionType`
  */
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
      internal::Region_growing_traits_2<Traits>,
      internal::Region_growing_traits_3<Traits> >::type;

    using Point = typename Polyline_traits::Point;
    using Vector = typename Polyline_traits::Vector;
    using Line = typename Polyline_traits::Line;

    using Squared_length = typename Polyline_traits::Compute_squared_length;
    using Squared_distance = typename Polyline_traits::Compute_squared_distance;
    using Scalar_product = typename Polyline_traits::Compute_scalar_product;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param input_range
      an instance of `InputRange` with polyline vertices

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \param point_map
      an instance of `PointMap` that maps an item from `input_range`
      to `Kernel::Point_2` or `Kernel::Point_3`

      \param traits
      an instance of `GeomTraits`

      \cgalNamedParamsBegin
        \cgalParamNBegin{distance_threshold}
          \cgalParamDescription{the maximum distance from a vertex to a line}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{1}
        \cgalParamNEnd
        \cgalParamNBegin{angle_threshold}
          \cgalParamDescription{the maximum accepted angle in degrees between
          the direction of a polyline edge and the direction of a line}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{25 degrees}
        \cgalParamNEnd
        \cgalParamNBegin{cos_value_threshold}
          \cgalParamDescription{the cos value computed as `cos(angle_threshold * PI / 180)`,
          this parameter can be used instead of the `angle_threshold`}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{`cos(25 * PI / 180)`}
        \cgalParamNEnd
        \cgalParamNBegin{min_region_size}
          \cgalParamDescription{the minimum number of vertices a region must have}
          \cgalParamType{`std::size_t`}
          \cgalParamDefault{2}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() > 0`
      \pre `distance_threshold >= 0`
      \pre `angle_threshold >= 0 && angle_threshold <= 90`
      \pre `cos_value_threshold >= 0 && cos_value_threshold <= 1`
      \pre `min_region_size > 0`
    */
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

      const FT angle_threshold = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::angle_threshold), FT(25));
      CGAL_precondition(angle_threshold >= FT(0) && angle_threshold <= FT(90));

      m_min_region_size = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::min_region_size), 2);
      CGAL_precondition(m_min_region_size > 0);

      const FT cos_value_threshold = static_cast<FT>(std::cos(CGAL::to_double(
        (angle_threshold * static_cast<FT>(CGAL_PI)) / FT(180))));
      m_cos_value_threshold = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::cos_value_threshold), cos_value_threshold);
      CGAL_precondition(m_cos_value_threshold >= FT(0) && m_cos_value_threshold <= FT(1));
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `RegionType::is_part_of_region()`.

      This function controls if a vertex with the index `query_index` is within
      the `distance_threshold` from the corresponding line and if the angle
      between the direction of the inward edge and the line's direction is within the `angle_threshold`.
      If both conditions are satisfied, it returns `true`, otherwise `false`.

      \param index1
      index of the previous vertex

      \param index2
      index of the query vertex

      \param region
      indices of vertices included in the region

      \return Boolean `true` or `false`

      \pre `region.size() > 0`
      \pre `index1 < input_range.size()`
      \pre `index2 < input_range.size()`
    */
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

      // Update new reference line and direction.
      if (m_direction_of_best_fit == CGAL::NULL_VECTOR) {
        if (input_point == query_point) return true;
        CGAL_precondition(input_point != query_point);
        m_line_of_best_fit = Line(input_point, query_point);
        m_direction_of_best_fit = m_line_of_best_fit.to_vector();
        return true;
      }
      CGAL_precondition(m_direction_of_best_fit != CGAL::NULL_VECTOR);

      // Add equal points to the previously defined region.
      if (input_point == query_point) return true;
      CGAL_precondition(input_point != query_point);
      const Vector query_direction(input_point, query_point);

      // Check real conditions.
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

    /*!
      \brief implements `RegionType::is_valid_region()`.

      This function controls if the `region` contains at least `min_region_size` vertices.

      \param region
      indices of vertices included in the region

      \return Boolean `true` or `false`
    */
    inline bool is_valid_region(const std::vector<std::size_t>& region) const {
      if (m_direction_of_best_fit == CGAL::NULL_VECTOR)
        return false; // all points are equal
      return (region.size() >= m_min_region_size);
    }

    /*!
      \brief implements `RegionType::update()`.

      This function fits the least squares line to all vertices from the `region`.

      \param region
      indices of vertices included in the region

      \return Boolean `true` if the line fitting succeeded and `false` otherwise

      \pre `region.size() > 0`
    */
    bool update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference line and direction
        m_direction_of_best_fit = CGAL::NULL_VECTOR;
      } else { // update reference line and direction
        if (m_direction_of_best_fit == CGAL::NULL_VECTOR)
          return false; // all points are equal
        CGAL_precondition(region.size() >= 2);
        std::tie(m_line_of_best_fit, m_direction_of_best_fit) =
          get_line_and_direction(region);
      }
      return true;
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    std::pair<Line, Vector> get_line_and_direction(
      const std::vector<std::size_t>& region) const {

      // The best fit line will be a line fitted to all region points with
      // its direction being the line's direction.
      CGAL_precondition(region.size() > 0);
      const Line line_of_best_fit =
        m_polyline_traits.create_line(
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
