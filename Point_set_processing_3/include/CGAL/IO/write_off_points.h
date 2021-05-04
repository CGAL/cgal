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

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Iterator_range.h>

#include <boost/utility/enable_if.hpp>

#include <iostream>
#include <fstream>
#include <fstream>
#include <iterator>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#endif

namespace CGAL {
namespace Point_set_processing_3 {
namespace internal {

template <typename PointRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_PSP(std::ostream& os,
                   const PointRange& points,
                   const CGAL_BGL_NP_CLASS& np)
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  // basic geometric types
  typedef typename CGAL::GetPointMap<PointRange, CGAL_BGL_NP_CLASS>::type                         PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, CGAL_BGL_NP_CLASS>::type      NormalMap;

  bool has_normals = !(std::is_same<NormalMap,
                                    typename Point_set_processing_3::GetNormalMap<
                                      PointRange, CGAL_BGL_NP_CLASS>::NoMap>::value);

  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  NormalMap normal_map = choose_parameter<NormalMap>(get_parameter(np, internal_np::normal_map));

  CGAL_point_set_processing_precondition(points.begin() != points.end());

  if(!os)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  set_stream_precision_from_NP(os, np);

  // Write header
  os << "NOFF" << std::endl;
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
      \cgalParamDefault{`the precision of the stream `os``}
    \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if writing was successful, `false` otherwise.
*/
template <typename PointRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(std::ostream& os,
               const PointRange& points,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::enable_if<internal::is_Range<PointRange> >::type* = nullptr
#endif
               )
{
  return Point_set_processing_3::internal::write_OFF_PSP(os, points, np);
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange>
bool write_OFF(std::ostream& os, const PointRange& points,
               typename boost::enable_if<internal::is_Range<PointRange> >::type* = nullptr)
{
  return write_OFF(os, points, parameters::all_default());
}

/// \endcond

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
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(const std::string& filename,
               const PointRange& points,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::enable_if<internal::is_Range<PointRange> >::type* = nullptr
#endif
               )
{
  std::ofstream os(filename);
  set_stream_precision_from_NP(os, np);
  return write_OFF(os, points, np);
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange>
bool write_OFF(const std::string& filename, const PointRange& points,
               typename boost::enable_if<internal::is_Range<PointRange> >::type* = nullptr)
{
  std::ofstream os(filename);
  return write_OFF(os, points, parameters::all_default());
}

/// \endcond

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

/**
  \ingroup PkgPointSetProcessing3IODeprecated

  \deprecated This function is deprecated since \cgal 5.3,
              \link PkgPointSetProcessing3IOOff `CGAL::IO::write_OFF()` \endlink should be used instead.
*/
template <typename PointRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool write_off_points(std::ostream& os, const PointRange& points, const CGAL_BGL_NP_CLASS& np)
{
  return IO::write_OFF(os, points, np);
}

/// \cond SKIP_IN_MANUAL

// variant with default NP
template <typename PointRange>
CGAL_DEPRECATED bool write_off_points(std::ostream& os, const PointRange& points)
{
  return IO::write_OFF(os, points, parameters::all_default());
}

/// \endcond

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_POINT_SET_PROCESSING_WRITE_OFF_POINTS_H
