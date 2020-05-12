// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_LINE_FIT_REGION_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_LINE_FIT_REGION_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>

namespace CGAL {
namespace Shape_detection {
namespace Point_set {

  /*!
    \ingroup PkgShapeDetectionRGOnPoints

    \brief Region type based on the quality of the least squares line
    fit applied to 2D points.

    This class fits a line, using \ref PkgPrincipalComponentAnalysisDRef "PCA",
    to chunks of points in a 2D point set and controls the quality of this fit.
    If all quality conditions are satisfied, the chunk is accepted as a valid region,
    otherwise rejected.

    \tparam GeomTraits
    must be a model of `Kernel`.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam PointMap
    must be an `LvaluePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Point_2`.

    \tparam NormalMap
    must be an `LvaluePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Vector_2`.

    \cgalModels `RegionType`
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap,
  typename NormalMap>
  class Least_squares_line_fit_region {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Normal_map = NormalMap;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// \cond SKIP_IN_MANUAL
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;

    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_FT = typename Local_traits::FT;
    using Local_point_2 = typename Local_traits::Point_2;
    using Local_line_2 = typename Local_traits::Line_2;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;

    using Squared_length_2 = typename Traits::Compute_squared_length_2;
    using Squared_distance_2 = typename Traits::Compute_squared_distance_2;
    using Scalar_product_2 = typename Traits::Compute_scalar_product_2;

    using Get_sqrt = internal::Get_sqrt<Traits>;
    using Sqrt = typename Get_sqrt::Sqrt;
    /// \endcond

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param input_range
      an instance of `InputRange` with 2D points and
      corresponding 2D normal vectors

      \param distance_threshold
      the maximum distance from a point to a line. %Default is 1.

      \param angle_threshold
      the maximum accepted angle in degrees between the normal of a point and
      the normal of a line. %Default is 25 degrees.

      \param min_region_size
      the minimum number of 2D points a region must have. %Default is 2.

      \param point_map
      an instance of `PointMap` that maps an item from `input_range`
      to `Kernel::Point_2`

      \param normal_map
      an instance of `NormalMap` that maps an item from `input_range`
      to `Kernel::Vector_2`

      \param traits
      an instance of `GeomTraits`

      \pre `input_range.size() > 0`
      \pre `distance_threshold >= 0`
      \pre `angle_threshold >= 0 && angle_threshold <= 90`
      \pre `min_region_size > 0`
    */
    Least_squares_line_fit_region(
      const InputRange& input_range,
      const FT distance_threshold = FT(1),
      const FT angle_threshold = FT(25),
      const std::size_t min_region_size = 2,
      const PointMap point_map = PointMap(),
      const NormalMap normal_map = NormalMap(),
      const GeomTraits traits = GeomTraits()) :
    m_input_range(input_range),
    m_distance_threshold(distance_threshold),
    m_normal_threshold(static_cast<FT>(
      std::cos(
        CGAL::to_double(
          (angle_threshold * static_cast<FT>(CGAL_PI)) / FT(180))))),
    m_min_region_size(min_region_size),
    m_point_map(point_map),
    m_normal_map(normal_map),
    m_squared_length_2(traits.compute_squared_length_2_object()),
    m_squared_distance_2(traits.compute_squared_distance_2_object()),
    m_scalar_product_2(traits.compute_scalar_product_2_object()),
    m_sqrt(Get_sqrt::sqrt_object(traits)),
    m_to_local_converter() {

      CGAL_precondition(input_range.size() > 0);

      CGAL_precondition(distance_threshold >= FT(0));
      CGAL_precondition(angle_threshold >= FT(0) && angle_threshold <= FT(90));
      CGAL_precondition(min_region_size > 0);
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `RegionType::is_part_of_region()`.

      This function controls if a point with the index `query_index` is within
      the `distance_threshold` from the corresponding line and if the angle
      between its normal and the line's normal is within the `angle_threshold`.
      If both conditions are satisfied, it returns `true`, otherwise `false`.

      \param query_index
      index of the query point

      The first and third parameters are not used in this implementation.

      \return Boolean `true` or `false`

      \pre `query_index >= 0 && query_index < input_range.size()`
    */
    bool is_part_of_region(
      const std::size_t,
      const std::size_t query_index,
      const std::vector<std::size_t>&) const {
      CGAL_precondition(query_index < m_input_range.size());

      const auto& key = *(m_input_range.begin() + query_index);
      const Point_2& query_point = get(m_point_map, key);

      const Vector_2& normal = get(m_normal_map, key);
      const FT normal_length = m_sqrt(m_squared_length_2(normal));
      CGAL_precondition(normal_length > FT(0));
      const Vector_2 query_normal = normal / normal_length;

      const FT distance_to_fitted_line =
      m_sqrt(m_squared_distance_2(query_point, m_line_of_best_fit));

      const FT cos_value =
      CGAL::abs(m_scalar_product_2(query_normal, m_normal_of_best_fit));

      return (( distance_to_fitted_line <= m_distance_threshold ) &&
        ( cos_value >= m_normal_threshold ));
    }

    /*!
      \brief implements `RegionType::is_valid_region()`.

      This function controls if the `region` contains at least `min_region_size` points.

      \param region
      indices of points included in the region

      \return Boolean `true` or `false`
    */
    inline bool is_valid_region(const std::vector<std::size_t>& region) const {
      return ( region.size() >= m_min_region_size );
    }

    /*!
      \brief implements `RegionType::update()`.

      This function fits the least squares line to all points from the `region`.

      \param region
      indices of points included in the region

      \pre `region.size() > 0`
    */
    void update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference line and normal
        CGAL_precondition(region[0] < m_input_range.size());

        // The best fit line will be a line through this point with
        // its normal being the point's normal.
        const auto& key = *(m_input_range.begin() + region[0]);

        const Point_2& point = get(m_point_map, key);
        const Vector_2& normal = get(m_normal_map, key);

        const FT normal_length = m_sqrt(m_squared_length_2(normal));

        CGAL_precondition(normal_length > FT(0));
        m_normal_of_best_fit =
        normal / normal_length;

        m_line_of_best_fit =
        Line_2(point, m_normal_of_best_fit).perpendicular(point);

      } else { // update reference line and normal

        std::vector<Local_point_2> points;
        points.reserve(region.size());

        for (std::size_t i = 0; i < region.size(); ++i) {
          CGAL_precondition(region[i] < m_input_range.size());

          const auto& key = *(m_input_range.begin() + region[i]);
          points.push_back(m_to_local_converter(get(m_point_map, key)));
        }
        CGAL_postcondition(points.size() == region.size());

        Local_line_2 fitted_line;
        Local_point_2 fitted_centroid;

        // The best fit line will be a line fitted to all region points with
        // its normal being perpendicular to the line.
        CGAL::linear_least_squares_fitting_2(
          points.begin(), points.end(),
          fitted_line, fitted_centroid,
          CGAL::Dimension_tag<0>(),
          Local_traits(),
          CGAL::Eigen_diagonalize_traits<Local_FT, 2>());

        m_line_of_best_fit =
        Line_2(
          static_cast<FT>(fitted_line.a()),
          static_cast<FT>(fitted_line.b()),
          static_cast<FT>(fitted_line.c()));

        const Vector_2 normal =
        m_line_of_best_fit.perpendicular(m_line_of_best_fit.point(0)).to_vector();
        const FT normal_length = m_sqrt(m_squared_length_2(normal));

        CGAL_precondition(normal_length > FT(0));
        m_normal_of_best_fit = normal / normal_length;
      }
    }

    /// @}

  private:

    // Fields.
    const Input_range& m_input_range;

    const FT m_distance_threshold;
    const FT m_normal_threshold;
    const std::size_t m_min_region_size;

    const Point_map m_point_map;
    const Normal_map m_normal_map;

    const Squared_length_2 m_squared_length_2;
    const Squared_distance_2 m_squared_distance_2;
    const Scalar_product_2 m_scalar_product_2;
    const Sqrt m_sqrt;

    const To_local_converter m_to_local_converter;

    Line_2 m_line_of_best_fit;
    Vector_2 m_normal_of_best_fit;
  };

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_LINE_FIT_REGION_H
