// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_WRITE_OFF_POINTS_H
#define CGAL_WRITE_OFF_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>


#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Iterator_range.h>

#include <iostream>
#include <iterator>

namespace CGAL {


/**
   \ingroup PkgPointSetProcessing3IO
   Saves the range of `points` (positions + normals, if available) to a .off ASCII stream.
   The function writes for each point a line with the x y z position
   followed by the nx ny nz normal (if available).

   \note The <A HREF="https://en.cppreference.com/w/cpp/io/ios_base/precision">`precision()`</A>
         of the output stream might not be sufficient depending on the data to be written.

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param stream output stream.
   \param points input point range.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the poing range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals are not written in the output stream.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \return `true` on success.
*/
template <typename PointRange,
          typename NamedParameters
>
bool
write_off_points(
  std::ostream& stream,
  const PointRange& points,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // basic geometric types
  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;

  bool has_normals = !(boost::is_same<NormalMap,
                       typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value);

  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  NormalMap normal_map = choose_parameter<NormalMap>(get_parameter(np, internal_np::normal_map));

  CGAL_point_set_processing_precondition(points.begin() != points.end());

  if(!stream)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // Write header
  stream << "NOFF" << std::endl;
  stream << points.size() << " 0 0" << std::endl;

  // Write positions + normals
  for(typename PointRange::const_iterator it = points.begin(); it != points.end(); it++)
  {
    stream << get(point_map, *it);
    if (has_normals)
      stream << " " << get(normal_map, *it);
    stream << std::endl;
  }

  return ! stream.fail();
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename PointRange>
bool
write_off_points(
  std::ostream& stream, ///< output stream.
  const PointRange& points)
{
  return write_off_points
    (stream, points, CGAL::Point_set_processing_3::parameters::all_default(points));
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template <typename ForwardIterator,
          typename PointMap,
          typename NormalMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points_and_normals(), please update your code")
bool
write_off_points_and_normals(
  std::ostream& stream, ///< output stream.
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3.
  NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_off_points
    (stream, points,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map).
     geom_traits(Kernel()));
}

// deprecated API
template <typename ForwardIterator,
          typename PointMap,
          typename NormalMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points_and_normals(), please update your code")
bool
write_off_points_and_normals(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond, ///< past-the-end input point.
  PointMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
  NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_off_points
    (stream, points,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map));
}

// deprecated API
template <typename ForwardIterator,
          typename NormalMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points_and_normals(), please update your code")
bool
write_off_points_and_normals(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond, ///< past-the-end input point.
  NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_off_points
    (stream, points,
     CGAL::parameters::normal_map (normal_map));
}

// deprecated API
template <typename ForwardIterator,
          typename PointMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points(), please update your code")
bool
write_off_points(
  std::ostream& stream, ///< output stream.
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3.
  const Kernel& ) ///< geometric traits.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_off_points
    (stream, points,
     CGAL::parameters::point_map (point_map).
     geom_traits (Kernel()));
}

// deprecated API
template <typename ForwardIterator,
          typename PointMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points(), please update your code")
bool
write_off_points(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond, ///< past-the-end input point.
  PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_off_points
    (stream, points,
     CGAL::parameters::point_map (point_map));
}

// deprecated API
template <typename ForwardIterator
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points(), please update your code")
bool
write_off_points(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond) ///< past-the-end input point.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_off_points
    (stream, points);
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond


} //namespace CGAL

#endif // CGAL_WRITE_OFF_POINTS_H
