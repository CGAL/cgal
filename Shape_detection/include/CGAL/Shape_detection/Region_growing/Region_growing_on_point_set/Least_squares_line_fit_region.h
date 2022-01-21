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
    a model of `Kernel`

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Point_2`

    \tparam NormalMap
    a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Vector_2`

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

    /// @}

  private:
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;

    using Squared_length_2 = typename Traits::Compute_squared_length_2;
    using Squared_distance_2 = typename Traits::Compute_squared_distance_2;
    using Scalar_product_2 = typename Traits::Compute_scalar_product_2;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param input_range
      an instance of `InputRange` with 2D points and
      corresponding 2D normal vectors

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{maximum_distance}
          \cgalParamDescription{the maximum distance from a point to a line}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{1}
        \cgalParamNEnd
        \cgalParamNBegin{maximum_angle}
          \cgalParamDescription{the maximum angle in degrees between
          the normal of a point and the normal of a line}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{25 degrees}
        \cgalParamNEnd
        \cgalParamNBegin{cosine_value}
          \cgalParamDescription{the cos value computed as `cos(maximum_angle * PI / 180)`,
          this parameter can be used instead of the `maximum_angle`}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{`cos(25 * PI / 180)`}
        \cgalParamNEnd
        \cgalParamNBegin{minimum_region_size}
          \cgalParamDescription{the minimum number of 2D points a region must have}
          \cgalParamType{`std::size_t`}
          \cgalParamDefault{2}
        \cgalParamNEnd
        \cgalParamNBegin{point_map}
          \cgalParamDescription{an instance of `PointMap` that maps an item from `input_range`
          to `Kernel::Point_2`}
          \cgalParamDefault{`PointMap()`}
        \cgalParamNEnd
        \cgalParamNBegin{normal_map}
          \cgalParamDescription{an instance of `NormalMap` that maps an item from `input_range`
          to `Kernel::Vector_2`}
          \cgalParamDefault{`NormalMap()`}
        \cgalParamNEnd
        \cgalParamNBegin{geom_traits}
          \cgalParamDescription{an instance of `GeomTraits`}
          \cgalParamDefault{`GeomTraits()`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() > 0`
      \pre `maximum_distance >= 0`
      \pre `maximum_angle >= 0 && maximum_angle <= 90`
      \pre `cosine_value >= 0 && cosine_value <= 1`
      \pre `minimum_region_size > 0`
    */
    template<typename NamedParameters>
    Least_squares_line_fit_region(
      const InputRange& input_range,
      const NamedParameters& np) :
    m_input_range(input_range),
    m_point_map(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::point_map), PointMap())),
    m_normal_map(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::normal_map), NormalMap())),
    m_traits(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::geom_traits), GeomTraits())),
    m_squared_length_2(m_traits.compute_squared_length_2_object()),
    m_squared_distance_2(m_traits.compute_squared_distance_2_object()),
    m_scalar_product_2(m_traits.compute_scalar_product_2_object()) {

      CGAL_precondition(input_range.size() > 0);
      const FT max_distance = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_distance), FT(1));
      CGAL_precondition(max_distance >= FT(0));
      m_distance_threshold = max_distance;

      const FT max_angle = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_angle), FT(25));
      CGAL_precondition(max_angle >= FT(0) && max_angle <= FT(90));

      m_min_region_size = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::minimum_region_size), 2);
      CGAL_precondition(m_min_region_size > 0);

      const FT default_cos_value = static_cast<FT>(std::cos(CGAL::to_double(
        (max_angle * static_cast<FT>(CGAL_PI)) / FT(180))));
      const FT cos_value = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::cosine_value), default_cos_value);
      CGAL_precondition(cos_value >= FT(0) && cos_value <= FT(1));
      m_cos_value_threshold = cos_value;
    }

    #if !defined(CGAL_NO_DEPRECATED_CODE) || defined(DOXYGEN_RUNNING)

    /*!
      \brief initializes all internal data structures.

      \deprecated This constructor is deprecated since the version 5.4 of \cgal.

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
    CGAL_DEPRECATED_MSG("This constructor is deprecated since the version 5.4 of CGAL!")
    Least_squares_line_fit_region(
      const InputRange& input_range,
      const FT distance_threshold = FT(1),
      const FT angle_threshold = FT(25),
      const std::size_t min_region_size = 2,
      const PointMap point_map = PointMap(),
      const NormalMap normal_map = NormalMap(),
      const GeomTraits traits = GeomTraits()) :
    Least_squares_line_fit_region(
      input_range, CGAL::parameters::
    maximum_distance(distance_threshold).
    maximum_angle(angle_threshold).
    minimum_region_size(min_region_size).
    point_map(point_map).
    normal_map(normal_map).
    geom_traits(traits))
    { }

    #endif // CGAL_NO_DEPRECATED_CODE

    /// \cond SKIP_IN_MANUAL
    // TODO: Should be off until the deprecated code is removed.
    // Least_squares_line_fit_region(
    //   const InputRange& input_range) :
    // Least_squares_line_fit_region(
    //   input_range, CGAL::parameters::all_default())
    // { }
    /// \endcond

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `RegionType::is_part_of_region()`.

      This function controls if a point with the index `query_index` is within
      the `maximum_distance` from the corresponding line and if the angle between
      its normal and the line's normal is within the `maximum_angle`. If both conditions
      are satisfied, it returns `true`, otherwise `false`.

      \param query_index
      index of the query point

      The first and third parameters are not used in this implementation.

      \return Boolean `true` or `false`

      \pre `query_index < input_range.size()`
    */
    bool is_part_of_region(
      const std::size_t,
      const std::size_t query_index,
      const std::vector<std::size_t>&) const {

      CGAL_precondition(query_index < m_input_range.size());
      const auto& key = *(m_input_range.begin() + query_index);
      const Point_2& query_point = get(m_point_map, key);
      const Vector_2& query_normal = get(m_normal_map, key);

      const FT a = CGAL::abs(m_line_of_best_fit.a());
      const FT b = CGAL::abs(m_line_of_best_fit.b());
      const FT c = CGAL::abs(m_line_of_best_fit.c());
      if (a == FT(0) && b == FT(0) && c == FT(0))
        return false;

      const FT squared_distance_to_fitted_line =
        m_squared_distance_2(query_point, m_line_of_best_fit);
      const FT squared_distance_threshold =
        m_distance_threshold * m_distance_threshold;

      const FT cos_value =
        m_scalar_product_2(query_normal, m_normal_of_best_fit);
      const FT squared_cos_value = cos_value * cos_value;

      FT squared_cos_value_threshold =
        m_cos_value_threshold * m_cos_value_threshold;
      squared_cos_value_threshold *= m_squared_length_2(query_normal);
      squared_cos_value_threshold *= m_squared_length_2(m_normal_of_best_fit);

      return (
        ( squared_distance_to_fitted_line <= squared_distance_threshold ) &&
        ( squared_cos_value >= squared_cos_value_threshold ));
    }

    /*!
      \brief implements `RegionType::is_valid_region()`.

      This function controls if the `region` contains at least `minimum_region_size` points.

      \param region
      indices of points included in the region

      \return Boolean `true` or `false`
    */
    inline bool is_valid_region(const std::vector<std::size_t>& region) const {
      return (region.size() >= m_min_region_size);
    }

    /*!
      \brief implements `RegionType::update()`.

      This function fits the least squares line to all points from the `region`.

      \param region
      indices of points included in the region

      \return Boolean `true` if the line fitting succeeded and `false` otherwise

      \pre `region.size() > 0`
    */
    bool update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference line and normal
        const std::size_t point_index = region[0];
        CGAL_precondition(point_index < m_input_range.size());

        // The best fit line will be a line through this point with
        // its normal being the point's normal.
        const auto& key = *(m_input_range.begin() + point_index);
        const Point_2& point = get(m_point_map, key);
        const Vector_2& normal = get(m_normal_map, key);
        if (normal == CGAL::NULL_VECTOR) return false;

        CGAL_precondition(normal != CGAL::NULL_VECTOR);
        m_line_of_best_fit = Line_2(point, normal).perpendicular(point);
        m_normal_of_best_fit = m_line_of_best_fit.perpendicular(
          m_line_of_best_fit.point(0)).to_vector();

      } else { // update reference line and normal
        CGAL_precondition(region.size() >= 2);
        std::tie(m_line_of_best_fit, m_normal_of_best_fit) =
          get_line_and_normal(region);
      }
      return true;
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    std::pair<Line_2, Vector_2> get_line_and_normal(
      const std::vector<std::size_t>& region) const {

      // The best fit line will be a line fitted to all region points with
      // its normal being perpendicular to the line.
      CGAL_precondition(region.size() > 0);
      const Line_2 unoriented_line_of_best_fit =
        internal::create_line_2(
          m_input_range, m_point_map, region, m_traits).first;
      const Vector_2 unoriented_normal_of_best_fit =
        unoriented_line_of_best_fit.perpendicular(
          unoriented_line_of_best_fit.point(0)).to_vector();

      // Flip the line's normal to agree with all input normals.
      long votes_to_keep_normal = 0;
      for (const std::size_t normal_index : region) {
        CGAL_precondition(normal_index < m_input_range.size());
        const auto& key = *(m_input_range.begin() + normal_index);
        const Vector_2& normal = get(m_normal_map, key);
        const bool agrees =
          m_scalar_product_2(normal, unoriented_normal_of_best_fit) > FT(0);
        votes_to_keep_normal += (agrees ? 1 : -1);
      }
      const bool flip_normal = (votes_to_keep_normal < 0);

      const Line_2 line_of_best_fit = flip_normal
        ? unoriented_line_of_best_fit.opposite()
        : unoriented_line_of_best_fit;
      const Vector_2 normal_of_best_fit = flip_normal
        ? (-1 * unoriented_normal_of_best_fit)
        : unoriented_normal_of_best_fit;

      return std::make_pair(line_of_best_fit, normal_of_best_fit);
    }
    /// \endcond

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const Normal_map m_normal_map;
    const Traits m_traits;

    FT m_distance_threshold;
    FT m_cos_value_threshold;
    std::size_t m_min_region_size;

    const Squared_length_2 m_squared_length_2;
    const Squared_distance_2 m_squared_distance_2;
    const Scalar_product_2 m_scalar_product_2;

    Line_2 m_line_of_best_fit;
    Vector_2 m_normal_of_best_fit;
  };

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_LINE_FIT_REGION_H
