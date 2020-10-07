// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_SCANLINE_ORIENT_NORMALS_H
#define CGAL_SCANLINE_ORIENT_NORMALS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/squared_distance_3.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <queue>

namespace CGAL
{

/// \cond SKIP_IN_MANUAL
namespace Point_set_processing_3
{

namespace internal
{

template <typename Iterator, typename PointMap,
          typename ScanDirectionFlagMap>
bool is_end_of_scanline (Iterator scanline_begin, Iterator it,
                         PointMap,
                         ScanDirectionFlagMap scan_direction_flag_map,
                         const Tag_false&) // no fallback
{
  return (get (scan_direction_flag_map, *scanline_begin)
          != get (scan_direction_flag_map, *it));
}

template <typename Iterator, typename PointMap,
          typename ScanDirectionFlagMap>
bool is_end_of_scanline (Iterator scanline_begin, Iterator it,
                         PointMap point_map,
                         ScanDirectionFlagMap,
                         const Tag_true&) // fallback
{
  using Point_3 = typename boost::property_traits<PointMap>::value_type;
  using Vector_3 = typename Kernel_traits<Point_3>::Kernel::Vector_3;

  if (std::distance (scanline_begin, it) < 3)
    return false;

  Iterator n_minus_1 = it; -- n_minus_1;
  Iterator n_minus_2 = n_minus_1; -- n_minus_2;
  Iterator n_minus_3 = n_minus_2; -- n_minus_3;

  const Point_3& p_minus_1 = get(point_map, *n_minus_1);
  const Point_3& p_minus_2 = get(point_map, *n_minus_2);
  const Point_3& p_minus_3 = get(point_map, *n_minus_3);

  // End of scanline reached if inversion of direction
  Vector_3 v32 (p_minus_3, p_minus_2);
  v32 = Vector_3 (v32.x(), v32.y(), 0);
  Vector_3 v21 (p_minus_2, p_minus_1);
  v21 = Vector_3 (v21.x(), v21.y(), 0);

  return (v32 * v21 < 0);
}

template <typename Iterator, typename PointMap>
typename Kernel_traits<typename boost::property_traits
                       <PointMap>::value_type>::Kernel::Vector_3
scanline_direction (Iterator begin, Iterator end,
                    PointMap point_map)
{
  using Vector_3 = typename Kernel_traits<typename boost::property_traits
                                          <PointMap>::value_type>::Kernel::Vector_3;
  Iterator last = end; -- last;
  Vector_3 direction (get (point_map, *begin),
                      get (point_map, *last));
  direction = Vector_3 (direction.x(), direction.y(), 0);
  direction = direction / CGAL::approximate_sqrt (direction.squared_length());
  return direction;
}

template <typename Iterator, typename PointMap, typename NormalMap,
          typename ScanDirectionFlagMap, typename ScanAngleMap>
void orient_scanline (Iterator begin, Iterator end,
                      PointMap point_map,
                      NormalMap normal_map,
                      ScanDirectionFlagMap scan_direction_flag,
                      ScanAngleMap scan_angle_map,
                      const Tag_false&, // no fallback direction flag
                      const Tag_false&) // no fallback scan angle
{
  using Vector_3 = typename boost::property_traits<NormalMap>::value_type;

  static const Vector_3 vertical (0, 0, 1);

  Vector_3 direction
    = Point_set_processing_3::internal::scanline_direction
    (begin, end, point_map);
  if (get (scan_direction_flag, *begin) == 1)
    direction = -direction;

  for (Iterator it = begin; it != end; ++ it)
  {
    double angle = CGAL_PI * double(get (scan_angle_map, *it)) / 180.;
    Vector_3 line_of_sight
      = direction * std::sin(angle) + vertical * std::cos(angle);
    const Vector_3& normal = get (normal_map, *it);
    if (line_of_sight * normal < 0)
      put (normal_map, *it, -normal);
  }
}

template <typename Iterator, typename PointMap, typename NormalMap,
          typename ScanDirectionFlagMap, typename ScanAngleMap>
void orient_scanline (Iterator begin, Iterator end,
                      PointMap point_map,
                      NormalMap normal_map,
                      ScanDirectionFlagMap,
                      ScanAngleMap scan_angle_map,
                      const Tag_true&, // no fallback direction flag
                      const Tag_false&) // fallback scan angle
{
  using Point_3 = typename boost::property_traits<PointMap>::value_type;
  using Vector_3 = typename boost::property_traits<NormalMap>::value_type;

  // Estimate scanner position:
  // above point with minimum scan angle
  float min_scan_angle = (std::numeric_limits<float>::max)();
  Iterator chosen = begin;
  double min_z = (std::numeric_limits<double>::max)();
  double max_z = -(std::numeric_limits<double>::max)();

  for (Iterator it = begin; it != end; ++ it)
  {
    const Point_3& p = get (point_map, *it);
    min_z = (std::min(min_z, p.z()));
    max_z = (std::max(max_z, p.z()));
    float scan_angle = get (scan_angle_map, *it);
    if (scan_angle < min_scan_angle)
    {
      min_scan_angle = scan_angle;
      chosen = it;
    }
  }

  const Point_3& chosen_point = get (point_map, *chosen);
  Point_3 scan_position (chosen_point.x(),
                         chosen_point.y(),
                         min_z + 10 * (max_z - min_z));

  for (Iterator it = begin; it != end; ++ it)
  {
    Vector_3 line_of_sight (get(point_map, *it), scan_position);
    const Vector_3& normal = get (normal_map, *it);
    if (line_of_sight * normal < 0)
      put (normal_map, *it, -normal);
  }
}

template <typename Iterator, typename PointMap, typename NormalMap,
          typename ScanDirectionFlagMap, typename ScanAngleMap,
          typename FallbackFlag>
void orient_scanline (Iterator begin, Iterator end,
                      PointMap point_map,
                      NormalMap normal_map,
                      ScanDirectionFlagMap,
                      ScanAngleMap,
                      const FallbackFlag&, // either fallback direction flag or not
                      const Tag_true&) // fallback scan angle
{
  using Point_3 = typename boost::property_traits<PointMap>::value_type;
  using Vector_3 = typename boost::property_traits<NormalMap>::value_type;

  // Estimate scanner position:
  // average XY-projected point, located above
  double mean_x = 0.;
  double mean_y = 0.;
  double min_z = (std::numeric_limits<double>::max)();
  double max_z = -(std::numeric_limits<double>::max)();
  std::size_t nb = 0;

  for (Iterator it = begin; it != end; ++ it)
  {
    const Point_3& p = get (point_map, *it);
    mean_x += p.x();
    mean_y += p.y();
    min_z = (std::min(min_z, p.z()));
    max_z = (std::max(max_z, p.z()));
    ++ nb;
  }

  Point_3 scan_position (mean_x / nb, mean_y / nb,
                         min_z + 10 * (max_z - min_z));

  for (Iterator it = begin; it != end; ++ it)
  {
    Vector_3 line_of_sight (get(point_map, *it), scan_position);
    const Vector_3& normal = get (normal_map, *it);
    if (line_of_sight * normal < 0)
      put (normal_map, *it, -normal);
  }
}

} // namespace internal

} // namespace Point_set_processing_3
/// \endcond

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessing3Algorithms

   Orients the normals of the range of `points` by estimating a line
   of sight and checking its consistency with the current normal orientation.

   \warning this function requires the input `points` to be ordered
   along scanlines aligned on the XY-plane. It is typically designed
   for 2.5D urban datasets acquired through, for example, airborne
   LIDAR devices.

   First, scanlines are estimated as subranges of `points` by
   iterating on `points`:

   - if the named parameter `scan_direction_flag` is provided, the
     range is cutted everytime the flag (which tells if the scanner
     was moving in the positive or negative direction) changes.

   - if no direction flag map is provided, a fallback method simply
     cuts the range everytime 3 consecutive points form an acute angle
     on the projected XY-plane. This fallback method gives suboptimal
     results.

   Then, the line of sight (estimated vector between a point and the
   position of the scanner at its time of acquisition) is estimated:

   - if both `scan_direction_flag` and `scan_angle` are provided, the
     line of sight can be directly computed as a combination of the 2D
     vector of the projected scanline, the vertical vector and the
     scan angle.

   - if only `scan_angle` is provided, then for each scanline, the
     position of the scanner is estimated as being above the point of
     the scanline which has the minimum scan angle absolute
     value. This method is less optimal than the one using
     `scan_direction_flag`.

   - if none of these property maps are provided, then for each
     scanline, the position of the scanner is estimated as being above
     of the barycenter of the points of the scanline projected on the
     XY-plane. This fallback method gives suboptimal results.

   Once the line of sight is estimated for each point, the normals are
   oriented by checking if, for one each, the line of sight and the
   normal vector give a positive scalar product. If they don't, then
   the normal vector is inverted.

   \note this method gives optimal results when `scan_direction_flag`
   and `scan_angle` are provided. Correct results may still be
   produced in the absence of either one or both of these properties,
   as long as the point set is ordered in 2.5D scanlines.

   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point set `points`}
       \cgalParamType{a model of `WritablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Vector_3`}
     \cgalParamNEnd

     \cgalParamNBegin{scan_angle_map}
       \cgalParamDescription{a property map associating the angle of acquisition to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `float`}
     \cgalParamNEnd

     \cgalParamNBegin{scan_direction_flag}
       \cgalParamDescription{a property map associating a flag describing the direction of the acquisition to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `unsigned char`}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd
*/
template <typename PointRange, typename NamedParameters>
void scanline_orient_normals (PointRange& points, const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using Iterator = typename PointRange::iterator;
  using PointMap = typename CGAL::GetPointMap<PointRange, NamedParameters>::type;
  using NormalMap = typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type;
  using Kernel = typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel;
  using Point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;
  using ScanAngleMap = typename Point_set_processing_3::GetScanAngleMap
    <PointRange, NamedParameters>::type;
  using ScanDirectionFlagMap = typename Point_set_processing_3::GetScanDirectionFlag
    <PointRange, NamedParameters>::type;

  CGAL_static_assertion_msg(!(std::is_same<NormalMap,
                              typename Point_set_processing_3::GetNormalMap
                              <PointRange, NamedParameters>::NoMap>::value),
                            "Error: no normal map");

  using Fallback_scan_angle
    = Boolean_tag
    <std::is_same<ScanAngleMap,
                  typename Point_set_processing_3::GetScanAngleMap
                  <PointRange, NamedParameters>::NoMap>::value>;
  using Fallback_scan_direction_flag
    = Boolean_tag
    <std::is_same<ScanDirectionFlagMap,
                  typename Point_set_processing_3::GetScanDirectionFlag
                  <PointRange, NamedParameters>::NoMap>::value>;


  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  NormalMap normal_map = choose_parameter<NormalMap>(get_parameter(np, internal_np::normal_map));
  ScanAngleMap scan_angle_map = choose_parameter<ScanAngleMap>
    (get_parameter(np, internal_np::scan_angle));
  ScanDirectionFlagMap scan_direction_flag_map = choose_parameter<ScanDirectionFlagMap>
    (get_parameter(np, internal_np::scan_direction_flag));

  std::size_t nb_scanlines = 1;

  Iterator scanline_begin = points.begin();
  for (Iterator it = points.begin(); it != points.end(); ++ it)
  {
    if (Point_set_processing_3::internal::is_end_of_scanline
        (scanline_begin, it, point_map, scan_direction_flag_map,
         Fallback_scan_direction_flag()))
    {
      Point_set_processing_3::internal::orient_scanline
        (scanline_begin, it, point_map, normal_map,
         scan_direction_flag_map, scan_angle_map,
         Fallback_scan_direction_flag(), Fallback_scan_angle());

      scanline_begin = it;
      ++ nb_scanlines;
    }
  }

  Point_set_processing_3::internal::orient_scanline
    (scanline_begin, points.end(), point_map, normal_map,
     scan_direction_flag_map, scan_angle_map,
     Fallback_scan_direction_flag(), Fallback_scan_angle());

  std::cerr << nb_scanlines << " scanline(s) identified (mean length = "
            << std::size_t(points.size() / double(nb_scanlines))
            << " point(s))" << std::endl;
}

template <typename PointRange>
void scanline_orient_normals (PointRange& points)
{
  return scanline_orient_normals (points,
                               CGAL::Point_set_processing_3::parameters::all_default(points));
}

} // namespace CGAL


#endif // CGAL_SCANLINE_ORIENT_NORMALS_H
