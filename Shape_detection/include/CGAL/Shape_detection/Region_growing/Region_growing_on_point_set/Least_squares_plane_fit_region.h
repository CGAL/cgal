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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_REGION_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_REGION_H

#include <CGAL/license/Shape_detection.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>

namespace CGAL {
namespace Shape_detection {
namespace Point_set {

  /*!
    \ingroup PkgShapeDetectionRGOnPoints

    \brief Region type based on the quality of the least squares plane
    fit applied to 3D points.

    This class fits a plane, using \ref PkgPrincipalComponentAnalysisDRef "PCA",
    to chunks of points in a 3D point set and controls the quality of this fit.
    If all quality conditions are satisfied, the chunk is accepted as a valid region,
    otherwise rejected.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Point_3`

    \tparam NormalMap
    a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Vector_3`

    \cgalModels `RegionType`
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap,
  typename NormalMap>
  class Least_squares_plane_fit_region {

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
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;

    using Squared_length_3 = typename Traits::Compute_squared_length_3;
    using Squared_distance_3 = typename Traits::Compute_squared_distance_3;
    using Scalar_product_3 = typename Traits::Compute_scalar_product_3;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param input_range
      an instance of `InputRange` with 3D points and
      corresponding 3D normal vectors

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{maximum_distance}
          \cgalParamDescription{the maximum distance from a point to a plane}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{1}
        \cgalParamNEnd
        \cgalParamNBegin{maximum_angle}
          \cgalParamDescription{the maximum angle in degrees between
          the normal of a point and the normal of a plane}
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
          \cgalParamDescription{the minimum number of 3D points a region must have}
          \cgalParamType{`std::size_t`}
          \cgalParamDefault{3}
        \cgalParamNEnd
        \cgalParamNBegin{point_map}
          \cgalParamDescription{an instance of `PointMap` that maps an item from `input_range`
          to `Kernel::Point_3`}
          \cgalParamDefault{`PointMap()`}
        \cgalParamNEnd
        \cgalParamNBegin{normal_map}
          \cgalParamDescription{ an instance of `NormalMap` that maps an item from `input_range`
          to `Kernel::Vector_3`}
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
    template<typename CGAL_NP_TEMPLATE_PARAMETERS>
    Least_squares_plane_fit_region(
      const InputRange& input_range,
      const CGAL_NP_CLASS& np = parameters::default_values()) :
      m_input_range(input_range),
      m_point_map(Point_set_processing_3_np_helper<InputRange, CGAL_NP_CLASS>::get_const_point_map(input_range, np)),
      m_normal_map(Point_set_processing_3_np_helper<InputRange, CGAL_NP_CLASS>::get_normal_map(input_range, np)),
      m_traits(parameters::choose_parameter(parameters::get_parameter(
        np, internal_np::geom_traits), GeomTraits())),
      m_squared_length_3(m_traits.compute_squared_length_3_object()),
      m_squared_distance_3(m_traits.compute_squared_distance_3_object()),
      m_scalar_product_3(m_traits.compute_scalar_product_3_object()) {

      CGAL_precondition(input_range.size() > 0);
      const FT max_distance = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_distance), FT(1));
      CGAL_precondition(max_distance >= FT(0));
      m_distance_threshold = max_distance;

      const FT max_angle = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_angle), FT(25));
      CGAL_precondition(max_angle >= FT(0) && max_angle <= FT(90));

      m_min_region_size = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::minimum_region_size), 3);
      CGAL_precondition(m_min_region_size > 0);

      const FT default_cos_value = static_cast<FT>(std::cos(CGAL::to_double(
        (max_angle * static_cast<FT>(CGAL_PI)) / FT(180))));
      const FT cos_value = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::cosine_value), default_cos_value);
      CGAL_precondition(cos_value >= FT(0) && cos_value <= FT(1));
      m_cos_value_threshold = cos_value;
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `RegionType::is_part_of_region()`.

      This function controls if a point with the index `query_index` is within
      the `maximum_distance` from the corresponding plane and if the angle between
      its normal and the plane's normal is within the `maximum_angle`. If both conditions
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
      const Point_3& query_point = get(m_point_map, key);
      const Vector_3& query_normal = get(m_normal_map, key);

      const FT a = CGAL::abs(m_plane_of_best_fit.a());
      const FT b = CGAL::abs(m_plane_of_best_fit.b());
      const FT c = CGAL::abs(m_plane_of_best_fit.c());
      const FT d = CGAL::abs(m_plane_of_best_fit.d());
      if (a == FT(0) && b == FT(0) && c == FT(0) && d == FT(0))
        return false;

      const FT squared_distance_to_fitted_plane =
        m_squared_distance_3(query_point, m_plane_of_best_fit);
      const FT squared_distance_threshold =
        m_distance_threshold * m_distance_threshold;

      const FT cos_value =
        m_scalar_product_3(query_normal, m_normal_of_best_fit);
      const FT squared_cos_value = cos_value * cos_value;

      FT squared_cos_value_threshold =
        m_cos_value_threshold * m_cos_value_threshold;
      squared_cos_value_threshold *= m_squared_length_3(query_normal);
      squared_cos_value_threshold *= m_squared_length_3(m_normal_of_best_fit);

      return (
        ( squared_distance_to_fitted_plane <= squared_distance_threshold ) &&
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

      This function fits the least squares plane to all points from the `region`.

      \param region
      indices of points included in the region

      \return Boolean `true` if the plane fitting succeeded and `false` otherwise

      \pre `region.size() > 0`
    */
    bool update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference plane and normal
        const std::size_t point_index = region[0];
        CGAL_precondition(point_index < m_input_range.size());

        // The best fit plane will be a plane through this point with
        // its normal being the point's normal.
        const auto& key = *(m_input_range.begin() + point_index);
        const Point_3& point = get(m_point_map, key);
        const Vector_3& normal = get(m_normal_map, key);
        if (normal == CGAL::NULL_VECTOR) return false;

        CGAL_precondition(normal != CGAL::NULL_VECTOR);
        m_plane_of_best_fit = Plane_3(point, normal);
        m_normal_of_best_fit = m_plane_of_best_fit.orthogonal_vector();

      } else { // update reference plane and normal
        if (region.size() < 3) return false;
        CGAL_precondition(region.size() >= 3);
        std::tie(m_plane_of_best_fit, m_normal_of_best_fit) =
          get_plane_and_normal(region);
      }
      return true;
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    std::pair<Plane_3, Vector_3> get_plane_and_normal(
      const std::vector<std::size_t>& region) const {

      // The best fit plane will be a plane fitted to all region points with
      // its normal being perpendicular to the plane.
      CGAL_precondition(region.size() > 0);
      const Plane_3 unoriented_plane_of_best_fit =
        internal::create_plane(
          m_input_range, m_point_map, region, m_traits).first;
      const Vector_3 unoriented_normal_of_best_fit =
        unoriented_plane_of_best_fit.orthogonal_vector();

      // Flip the plane's normal to agree with all input normals.
      long votes_to_keep_normal = 0;
      for (const std::size_t normal_index : region) {
        CGAL_precondition(normal_index < m_input_range.size());
        const auto& key = *(m_input_range.begin() + normal_index);
        const Vector_3& normal = get(m_normal_map, key);
        const bool agrees =
          m_scalar_product_3(normal, unoriented_normal_of_best_fit) > FT(0);
        votes_to_keep_normal += (agrees ? 1 : -1);
      }
      const bool flip_normal = (votes_to_keep_normal < 0);

      const Plane_3 plane_of_best_fit = flip_normal
        ? unoriented_plane_of_best_fit.opposite()
        : unoriented_plane_of_best_fit;
      const Vector_3 normal_of_best_fit = flip_normal
        ? (-1 * unoriented_normal_of_best_fit)
        : unoriented_normal_of_best_fit;

      return std::make_pair(plane_of_best_fit, normal_of_best_fit);
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

    const Squared_length_3 m_squared_length_3;
    const Squared_distance_3 m_squared_distance_3;
    const Scalar_product_3 m_scalar_product_3;

    Plane_3 m_plane_of_best_fit;
    Vector_3 m_normal_of_best_fit;
  };

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_PLANE_FIT_REGION_H
