// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_CYLINDER_FIT_REGION_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_CYLINDER_FIT_REGION_H

#include <CGAL/license/Shape_detection.h>

#include <vector>

#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Eigen_diagonalize_traits.h>

#include <CGAL/Shape_detection/Region_growing/internal/utils.h>

namespace CGAL {
namespace Shape_detection {
namespace Point_set {

/*!
  \ingroup PkgShapeDetectionRGOnPoints

  \brief Region type based on the quality of the least squares cylinder
  fit applied to 3D points.

  This class fits an infinite cylinder to chunks of points in a 3D
  point set and controls the quality of this fit.  If all quality
  conditions are satisfied, the chunk is accepted as a valid region,
  otherwise rejected.

  \tparam GeomTraits
  must be a model of `Kernel`.

  \tparam InputRange
  must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

  \tparam PointMap
  must be an `LvaluePropertyMap` whose key type is the value type of the input
  range and value type is `Kernel::Point_3`.

  \tparam NormalMap
  must be an `LvaluePropertyMap` whose key type is the value type of the input
  range and value type is `Kernel::Vector_3`.

  \cgalModels `RegionType`
*/
template<typename GeomTraits,
         typename InputRange,
         typename PointMap,
         typename NormalMap>
class Least_squares_cylinder_fit_region
{

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
  using Point_3 = typename Traits::Point_3;
  using Vector_3 = typename Traits::Vector_3;
  using Line_3 = typename Traits::Line_3;

  using Squared_distance_3 = typename Traits::Compute_squared_distance_3;

  using Get_sqrt = internal::Get_sqrt<Traits>;
  using Sqrt = typename Get_sqrt::Sqrt;

private:

  const Input_range& m_input_range;

  const FT m_distance_threshold;
  const FT m_normal_threshold;
  const std::size_t m_min_region_size;
  const FT m_min_radius;
  const FT m_max_radius;

  const Point_map m_point_map;
  const Normal_map m_normal_map;

  const Squared_distance_3 m_squared_distance_3;
  const Sqrt m_sqrt;

  Line_3 m_axis;
  FT m_radius;

public:

  /// \endcond

  /// @}

  /// \name Initialization
  /// @{

  /*!
    \brief initializes all internal data structures.

    \param input_range an instance of `InputRange` with 3D points and
    corresponding 3D normal vectors

    \param distance_threshold the maximum distance from a point to a
    cylinder. %Default is 1.

    \param angle_threshold the maximum accepted angle in degrees
    between the normal of a point and the radius of a cylinder. %Default
    is 25 degrees.

    \param min_region_size the minimum number of 3D points a region
    must have. %Default is 3.

    \param minimum_radius is the radius below which an estimated
    cylinder is considered as invalid and discarded. %Default is 0 (no
    limit).

    \param maximum_radius is the radius above which an estimated
    cylinder is considered as invalid and discarded. %Default is
    infinity (no limit).

    \param point_map an instance of `PointMap` that maps an item from
    `input_range` to `Kernel::Point_3`

    \param normal_map an instance of `NormalMap` that maps an item
    from `input_range` to `Kernel::Vector_3`

    \param traits an instance of `GeomTraits`.

    \pre `input_range.size() > 0`
    \pre `distance_threshold >= 0`
    \pre `angle_threshold >= 0 && angle_threshold <= 90`
    \pre `min_region_size > 0`
  */
  Least_squares_cylinder_fit_region  (const InputRange& input_range,
                                      const FT distance_threshold = FT(1),
                                      const FT angle_threshold = FT(25),
                                      const std::size_t min_region_size = 3,
                                      const FT minimum_radius = FT(0),
                                      const FT maximum_radius = std::numeric_limits<FT>::infinity(),
                                      const PointMap point_map = PointMap(),
                                      const NormalMap normal_map = NormalMap(),
                                      const GeomTraits traits = GeomTraits())
  : m_input_range(input_range)
  , m_distance_threshold(distance_threshold)
  , m_normal_threshold(static_cast<FT>(
                         std::cos(
                           CGAL::to_double(
                             (angle_threshold * static_cast<FT>(CGAL_PI)) / FT(180)))))
  , m_min_region_size(min_region_size)
  , m_min_radius (minimum_radius)
  , m_max_radius (maximum_radius)
  , m_point_map(point_map)
  , m_normal_map(normal_map)
  , m_squared_distance_3(traits.compute_squared_distance_3_object())
  , m_sqrt(Get_sqrt::sqrt_object(traits))
  , m_radius(std::numeric_limits<FT>::quiet_NaN())
  {
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

    This function controls if a point with the index `query_index` is
    within the `distance_threshold` from the corresponding cylinder and
    if the angle between its normal and the cylinder radius is within
    the `angle_threshold`.  If both conditions are satisfied, it
    returns `true`, otherwise `false`.

    \param query_index
    index of the query point

    The first and third parameters are not used in this implementation.

    \return Boolean `true` or `false`

    \pre `query_index >= 0 && query_index < input_range.size()`
  */
  bool is_part_of_region (const std::size_t,
                          const std::size_t query_index,
                          const std::vector<std::size_t>& indices) const
  {
    CGAL_precondition(query_index < m_input_range.size());

    // First, we need to integrate at least 6 points so that the
    // computed cylinder means something
    if (indices.size() < 6)
      return true;

    if (std::isnan(m_radius))
      return false;

    // If radius is out of bound, nothing fits, early ending
    if (m_radius < m_min_radius || m_radius > m_max_radius)
      return false;

    const auto& key = *(m_input_range.begin() + query_index);
    const Point_3& query_point = get(m_point_map, key);
    Vector_3 normal = get(m_normal_map, key);

    FT distance_to_center = m_sqrt(m_squared_distance_3 (query_point, m_axis));
    FT distance_to_cylinder = CGAL::abs (distance_to_center - m_radius);

    if (distance_to_cylinder > m_distance_threshold)
      return false;

    normal = normal / m_sqrt (normal * normal);

    Point_3 proj = m_axis.projection(query_point);
    Vector_3 ray (proj, query_point);
    ray = ray / m_sqrt (ray * ray);

    if (CGAL::abs (normal * ray) < m_normal_threshold)
      return false;

    return true;
  }

  /*!
    \brief implements `RegionType::is_valid_region()`.

    This function controls if the estimated radius is between
    `minimum_radius` and `maximum_radius` and if the `region` contains
    at least `min_region_size` points.

    \param region
    indices of points included in the region

    \return Boolean `true` or `false`
  */
  inline bool is_valid_region(const std::vector<std::size_t>& region) const
  {
    return ((m_min_radius <= m_radius && m_radius <= m_max_radius)
            && (region.size() >= m_min_region_size));
  }

  /*!
    \brief implements `RegionType::update()`.

    This function fits the least squares cylinder to all points from the `region`.

    \param region
    indices of points included in the region

    \pre `region.size() > 0`
  */
  void update(const std::vector<std::size_t>& region)
  {
    CGAL_precondition(region.size() > 0);

    Vector_3 mean_axis = CGAL::NULL_VECTOR;
    std::size_t nb = 0;
    // Shuffle to avoid always picking 2 close points
    std::vector<std::size_t>& aregion
      = const_cast<std::vector<std::size_t>&>(region);
    cpp98::random_shuffle (aregion.begin(), aregion.end());

    const Point_3& ref = get(m_point_map, *(m_input_range.begin() + region.front()));

    m_radius = FT(0);
    Point_3 point_on_axis = CGAL::ORIGIN;
    for (std::size_t i = 0; i < aregion.size() - 1; ++ i)
    {
      Vector_3 v0 = get(m_normal_map, *(m_input_range.begin() + aregion[i]));
      v0 = v0 / CGAL::sqrt(v0*v0);
      Vector_3 v1 = get(m_normal_map, *(m_input_range.begin() + aregion[i+1]));
      v1 = v1 / CGAL::sqrt(v1*v1);
      Vector_3 axis = CGAL::cross_product (v0, v1);
      if (CGAL::sqrt(axis.squared_length()) < (FT)(0.01))
        continue;
      axis = axis / CGAL::sqrt(axis * axis);

      const Point_3& p0 = get(m_point_map, *(m_input_range.begin() + aregion[i]));
      const Point_3& p1 = get(m_point_map, *(m_input_range.begin() + aregion[i+1]));

      Vector_3 xdir = v0 - axis * (v0 * axis);
      xdir = xdir / CGAL::sqrt (xdir * xdir);

      Vector_3 ydir = CGAL::cross_product (axis, xdir);
      ydir = ydir / CGAL::sqrt (ydir * ydir);

      FT v1x = v1 * ydir;
      FT v1y = -v1 * xdir;

      Vector_3 d (p0, p1);

      FT ox = xdir * d;
      FT oy = ydir * d;
      FT ldist = v1x * ox + v1y * oy;

      FT radius = ldist / v1x;

      Point_3 point = p0 + xdir * radius;
      Line_3 line (point, axis);
      point = line.projection(ref);

      point_on_axis = CGAL::barycenter (point_on_axis, nb, point, 1);

      m_radius += CGAL::abs(radius);

      if (nb != 0 && axis * mean_axis < 0)
        axis = -axis;

      mean_axis = mean_axis + axis;
      ++ nb;
    }

    if (nb == 0)
      return;

    mean_axis = mean_axis / CGAL::sqrt(mean_axis * mean_axis);
    m_axis = Line_3 (point_on_axis, mean_axis);
    m_radius /= nb;

    m_radius = FT(0);
    for (std::size_t i = 0; i < aregion.size(); ++ i)
    {
      const Point_3& p0 = get(m_point_map, *(m_input_range.begin() + aregion[i]));
      m_radius += m_sqrt(m_squared_distance_3(p0, m_axis));
    }
    m_radius /= aregion.size();
  }
  /// @}

};

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_CYLINDER_FIT_REGION_H
