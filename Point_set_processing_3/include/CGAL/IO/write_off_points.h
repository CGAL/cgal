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

#ifndef CGAL_POINT_SET_PROCESSING_WRITE_OFF_POINTS_H
#define CGAL_POINT_SET_PROCESSING_WRITE_OFF_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/IO/helpers.h>
#include <CGAL/IO/OFF.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/assertions.h>

#include <iostream>
#include <fstream>
#include <fstream>
#include <iterator>
#include <type_traits>

namespace CGAL {
namespace Point_set_processing_3 {
namespace internal {

template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_PSP(std::ostream& os,
                   const PointRange& points,
                   const CGAL_NP_CLASS& np = CGAL::parameters::default_values())
{
  using CGAL::parameters::is_default_parameter;

  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, CGAL_NP_CLASS> NP_helper;
  typedef typename NP_helper::Const_point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;

  const bool has_normals = NP_helper::has_normal_map(points, np);

  PointMap point_map = NP_helper::get_const_point_map(points, np);
  NormalMap normal_map = NP_helper::get_normal_map(points, np);

  CGAL_precondition(points.begin() != points.end());

  if(!os)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  set_stream_precision_from_NP(os, np);

  // Write header
  if (has_normals)
    os << "NOFF" << std::endl;
  else
    os << "OFF" << std::endl;
  os << points.size() << " 0 0" << std::endl;

  // Write positions + normals
  for(typename PointRange::const_iterator it = points.begin(); it != points.end(); it++)
  {
    os << get(point_map, *it);
    if (has_normals)
      os << " " << get(normal_map, *it);
    os << "\n";
  }

  os << std::flush;

  return !os.fail();
}

} // namespace internal
} // namespace Point_set_processing_3

namespace IO {

/**
   \ingroup PkgPointSetProcessing3IOOff

   \brief writes the range of `points` (positions + normals, if available), using the \ref IOStreamOFF.

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param os output stream
   \param points input point range
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals are not written in the output stream.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd

    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{the precision of the stream `os`}
    \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if writing was successful, `false` otherwise.
*/
template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(std::ostream& os,
               const PointRange& points,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PointRange>::value>* = nullptr
#endif
               )
{
  return Point_set_processing_3::internal::write_OFF_PSP(os, points, np);
}

/**
   \ingroup PkgPointSetProcessing3IOOff

   \brief writes the range of `points` (positions + normals, if available), using the \ref IOStreamOFF.

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param filename the path to the output file
   \param points input point range
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals are not written in the output file.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd

    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`6`}
    \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if writing was successful, `false` otherwise.

   \sa \ref IOStreamOFF
*/
template <typename PointRange,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(const std::string& filename,
               const PointRange& points,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PointRange>::value>* = nullptr
#endif
               )
{
  std::ofstream os(filename);
  set_stream_precision_from_NP(os, np);
  return write_OFF(os, points, np);
}

} // IO namespace

#ifndef CGAL_NO_DEPRECATED_CODE

/// \cond SKIP_IN_MANUAL

template <typename ForwardIterator,
          typename PointMap,
          typename NormalMap,
          typename Kernel>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points_and_normals(), please update your code")
bool write_off_points_and_normals(std::ostream& os, ///< output stream.
                                  ForwardIterator first,  ///< iterator over the first input point.
                                  ForwardIterator beyond, ///< past-the-end iterator over the input points.
                                  PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3.
                                  NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3.
                                  const Kernel& /*kernel*/) ///< geometric traits.
{
  CGAL::Iterator_range<ForwardIterator> points(first, beyond);
  return write_off_points(os, points,
                          CGAL::parameters::point_map(point_map)
                                           .normal_map(normal_map)
                                           .geom_traits(Kernel()));
}

template <typename ForwardIterator,
          typename PointMap,
          typename NormalMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points_and_normals(), please update your code")
bool write_off_points_and_normals(std::ostream& os, ///< output stream.
                                  ForwardIterator first, ///< first input point.
                                  ForwardIterator beyond, ///< past-the-end input point.
                                  PointMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                  NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  CGAL::Iterator_range<ForwardIterator> points(first, beyond);
  return write_off_points(os, points,
                          parameters::point_map(point_map)
                                     .normal_map(normal_map));
}

template <typename ForwardIterator,
          typename NormalMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points_and_normals(), please update your code")
bool write_off_points_and_normals(std::ostream& os, ///< output stream.
                                  ForwardIterator first, ///< first input point.
                                  ForwardIterator beyond, ///< past-the-end input point.
                                  NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  CGAL::Iterator_range<ForwardIterator> points(first, beyond);
  return write_off_points(os, points, parameters::normal_map (normal_map));
}

template <typename ForwardIterator,
          typename PointMap,
          typename Kernel>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points(), please update your code")
bool write_off_points(std::ostream& os, ///< output stream.
                      ForwardIterator first,  ///< iterator over the first input point.
                      ForwardIterator beyond, ///< past-the-end iterator over the input points.
                      PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3.
                      const Kernel&) ///< geometric traits.
{
  CGAL::Iterator_range<ForwardIterator> points(first, beyond);
  return write_off_points(os, points,
                          parameters::point_map(point_map)
                                     .geom_traits(Kernel()));
}

template <typename ForwardIterator,
          typename PointMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points(), please update your code")
bool write_off_points(std::ostream& os, ///< output stream.
                      ForwardIterator first, ///< first input point.
                      ForwardIterator beyond, ///< past-the-end input point.
                      PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  CGAL::Iterator_range<ForwardIterator> points(first, beyond);
  return write_off_points(os, points, parameters::point_map (point_map));
}

template <typename ForwardIterator
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_off_points(), please update your code")
bool write_off_points(std::ostream& os, ///< output stream.
                      ForwardIterator first, ///< first input point.
                      ForwardIterator beyond) ///< past-the-end input point.
{
  CGAL::Iterator_range<ForwardIterator> points(first, beyond);
  return write_off_points(os, points);
}

/// \endcond

template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool write_off_points(std::ostream& os, const PointRange& points, const CGAL_NP_CLASS& np = parameters::default_values())
{
  return IO::write_OFF(os, points, np);
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_POINT_SET_PROCESSING_WRITE_OFF_POINTS_H
