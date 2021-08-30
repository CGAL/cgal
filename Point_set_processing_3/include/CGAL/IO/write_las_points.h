// Copyright (c) 2017  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_POINT_SET_PROCESSING_WRITE_LAS_POINTS_H
#define CGAL_POINT_SET_PROCESSING_WRITE_LAS_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/IO/helpers.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Kernel_traits.h>

#include <boost/cstdint.hpp>
#include <boost/version.hpp>
#include <boost/utility/enable_if.hpp>

#ifdef BOOST_MSVC
#  pragma warning(push)
#  pragma warning(disable:4251) // DLL warning from LASlib
#endif

#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

#define USE_AS_DLL
#include <lasdefinitions.hpp>
#include <lasreader_las.hpp>
#include <laswriter_las.hpp>
#undef USE_AS_DLL

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif

#ifdef BOOST_MSVC
#  pragma warning(pop)
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#define CGAL_DEPRECATED
#endif

namespace CGAL {

namespace IO {

/**
 \ingroup PkgPointSetProcessing3IOLas

 \brief generates a %LAS property handler to write 3D points.

 \tparam PointMap the property map used to store points.

 \sa `write_LAS()`
 \sa \ref IOStreamLAS
*/
template <typename PointMap>
std::tuple<PointMap, LAS_property::X, LAS_property::Y, LAS_property::Z >
make_las_point_writer(PointMap point_map)
{
  return std::make_tuple (point_map, LAS_property::X(), LAS_property::Y(), LAS_property::Z());
}

/// \cond SKIP_IN_MANUAL

namespace internal {
namespace LAS {

  inline void output_value(LASpoint& r, const unsigned short& v, LAS_property::Intensity&)
  { r.set_intensity(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, LAS_property::Return_number&)
  { r.set_return_number(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, LAS_property::Number_of_returns&)
  { r.set_number_of_returns(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, LAS_property::Scan_direction_flag&)
  { r.set_scan_direction_flag(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, LAS_property::Edge_of_flight_line&)
  { r.set_edge_of_flight_line(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, LAS_property::Classification&)
  { r.set_classification(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, LAS_property::Synthetic_flag&)
  { r.set_synthetic_flag(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, LAS_property::Keypoint_flag&)
  { r.set_keypoint_flag(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, LAS_property::Withheld_flag&)
  { r.set_withheld_flag(v); }
  inline void output_value(LASpoint& r, const float& v, LAS_property::Scan_angle&)
  { r.set_scan_angle_rank(char(v)); }
  inline void output_value(LASpoint& r, const unsigned char& v, LAS_property::User_data&)
  { r.set_user_data(v); }
  inline void output_value(LASpoint& r, const unsigned short& v, LAS_property::Point_source_ID&)
  { r.set_point_source_ID(v); }
  inline void output_value(LASpoint& r, const unsigned int& v, LAS_property::Deleted_flag&)
  { r.set_deleted_flag(v); }
  inline void output_value(LASpoint& r, const double& v, LAS_property::GPS_time&)
  { r.set_gps_time(v); }
  inline void output_value(LASpoint& r, const unsigned short& v, LAS_property::R&)
  { r.set_R(v); }
  inline void output_value(LASpoint& r, const unsigned short& v, LAS_property::G&)
  { r.set_G(v); }
  inline void output_value(LASpoint& r, const unsigned short& v, LAS_property::B&)
  { r.set_B(v); }
  inline void output_value(LASpoint& r, const unsigned short& v, LAS_property::I&)
  { r.set_I(v); }

  template <typename ForwardIterator>
  void output_properties(LASpoint&, ForwardIterator) { }

  template <typename ForwardIterator,
            typename PropertyMap,
            typename T>
  void output_properties(LASpoint& point,
                         ForwardIterator it,
                         std::pair<PropertyMap, T>&& current)
  {
    output_value (point, get(current.first, *it), current.second);
  }

  template <typename ForwardIterator,
            typename PropertyMap,
            typename T,
            typename NextPropertyHandler,
            typename ... PropertyHandler>
  void output_properties(LASpoint& point,
                         ForwardIterator it,
                         std::pair<PropertyMap, T>&& current,
                         NextPropertyHandler&& next,
                         PropertyHandler&& ... properties)
  {
    output_value (point, get(current.first, *it), current.second);
    output_properties (point, it, std::forward<NextPropertyHandler>(next),
                       std::forward<PropertyHandler>(properties)...);
  }

} // namespace LAS
} // namespace internal

/// \endcond

/**
   \ingroup PkgPointSetProcessing3IOLas

   \brief writes the range of `points` with properties to a .las stream.

   Properties are handled through a variadic list of property
   handlers. A `PropertyHandle` is a `std::pair<PropertyMap,
   LAS_property::Tag >` used to write a scalar value
   `LAS_property::Tag::type` as a %LAS property (for example,
   writing an `int` vairable as an `int` %LAS property). An exception
   is used for points that are written using a `std::tuple` object.

   See documentation of `read_LAS_with_properties()` for the
   list of available `LAS_property::Tag` classes.

   \attention When writing a binary file, the flag `std::ios::binary` flag must be set during the creation of the `ofstream`.

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam PointMap is a model of `ReadablePropertyMap` with a value_type = `CGAL::Point_3`.
   \tparam PropertyHandler handlers to recover properties.

   \returns `true` if writing was successful, `false` otherwise.

   \sa `make_las_point_writer()`
   \sa \ref IOStreamLAS
*/
template <typename PointRange,
          typename PointMap,
          typename ... PropertyHandler>
bool write_LAS_with_properties(std::ostream& os, ///< output stream.
                               const PointRange& points, ///< input point range.
                               std::tuple<PointMap,
                               LAS_property::X,
                               LAS_property::Y,
                               LAS_property::Z> point_property, ///< property handler for points
                               PropertyHandler&& ... properties) ///< parameter pack of property handlers
{
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  if(!os)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  CGAL::Bbox_3 bbox = CGAL::bbox_3
    (boost::make_transform_iterator
     (points.begin(), CGAL::Property_map_to_unary_function<PointMap>(std::get<0>(point_property))),
     boost::make_transform_iterator
     (points.end(), CGAL::Property_map_to_unary_function<PointMap>(std::get<0>(point_property))));

  LASheader header;
  header.x_scale_factor = 1e-9 * (bbox.xmax() - bbox.xmin());
  header.y_scale_factor = 1e-9 * (bbox.ymax() - bbox.ymin());
  header.z_scale_factor = 1e-9 * (bbox.zmax() - bbox.zmin());
  header.x_offset = bbox.xmin();
  header.y_offset = bbox.ymin();
  header.z_offset = bbox.zmin();
  header.point_data_format = 3;
  header.point_data_record_length = 34;

  LASpoint laspoint;
  laspoint.init(&header, header.point_data_format, header.point_data_record_length, 0);

  LASwriterLAS laswriter;
  laswriter.open (os, &header);

  // Write positions + normals
  for(typename PointRange::const_iterator it = points.begin(); it != points.end(); it++)
  {
    const typename PointMap::value_type& p = get(std::get<0>(point_property), *it);
    laspoint.set_X ((unsigned int)((p.x() - header.x_offset) / header.x_scale_factor));
    laspoint.set_Y ((unsigned int)((p.y() - header.y_offset) / header.y_scale_factor));
    laspoint.set_Z ((unsigned int)((p.z() - header.z_offset) / header.z_scale_factor));
    internal::LAS::output_properties (laspoint, it, std::forward<PropertyHandler>(properties)...);

    laswriter.write_point (&laspoint);
    laswriter.update_inventory(&laspoint);
  }

  laswriter.update_header(&header, TRUE);

  laswriter.close();

  return !os.fail();
}

/**
   \ingroup PkgPointSetProcessing3IOLas

   \brief writes the range of `points` (positions only), using the \ref IOStreamLAS.

   \attention When writing a binary file, the flag `std::ios::binary` flag must be set during the creation of the `ofstream`.

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

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if writing was successful, `false` otherwise.

   \sa \ref IOStreamLAS
   \sa `write_LAS_with_properties()`
*/
template <typename PointRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_LAS(std::ostream& os,
               const PointRange& points,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::enable_if<internal::is_Range<PointRange> >::type* = nullptr
#endif
               )
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetPointMap<PointRange, CGAL_BGL_NP_CLASS>::type PointMap;
  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));

  if(!os)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  return write_LAS_with_properties(os, points, make_las_point_writer(point_map));
}

/**
   \ingroup PkgPointSetProcessing3IOLas

   Saves the range of `points` (positions only), using the \ref IOStreamLAS.

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param filename the path the output file
   \param points input point range
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if writing was successful, `false` otherwise.

   \sa `write_LAS_with_properties()`
*/
template <typename PointRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_LAS(const std::string& filename,
               const PointRange& points,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::enable_if<internal::is_Range<PointRange> >::type* = nullptr
#endif
               )
{
  std::ofstream os(filename, std::ios::binary);
  CGAL::IO::set_mode(os, CGAL::IO::BINARY);
  return write_LAS(os, points, np);
}

/// \cond SKIP_IN_MANUAL

// variant with default NP
template <typename PointRange>
bool write_LAS(std::ostream& os, const PointRange& points,
               typename boost::enable_if<internal::is_Range<PointRange> >::type* = nullptr)
{
  return write_LAS(os, points, CGAL::Point_set_processing_3::parameters::all_default(points));
}

template <typename PointRange>
bool write_LAS(const std::string& filename, const PointRange& points,
               typename boost::enable_if<internal::is_Range<PointRange> >::type* = nullptr)
{
  std::ofstream os(filename, std::ios::binary);
  CGAL::IO::set_mode(os, CGAL::IO::BINARY);
  return write_LAS(os, points, parameters::all_default());
}

/// \endcond

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

using IO::make_las_point_writer;

/// \cond SKIP_IN_MANUAL

template <typename ForwardIterator,
          typename PointMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_las_points(), please update your code")
bool write_las_points(std::ostream& os, ///< output stream.
                      ForwardIterator first, ///< first input point.
                      ForwardIterator beyond, ///< past-the-end input point.
                      PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return IO::write_LAS(os, points, parameters::point_map(point_map));
}

template <typename ForwardIterator>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_las_points(), please update your code")
bool write_las_points(std::ostream& os, ///< output stream.
                      ForwardIterator first, ///< first input point.
                      ForwardIterator beyond) ///< past-the-end input point.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return IO::write_LAS(os, points);
}

/// \endcond

/**
  \ingroup PkgPointSetProcessing3IODeprecated

  \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::write_LAS_with_properties()` should be used instead.
*/
template <typename PointRange,
          typename PointMap,
          typename ... PropertyHandler>
CGAL_DEPRECATED bool write_las_points_with_properties(std::ostream& os,
                                                      const PointRange& points,
                                                      std::tuple<PointMap,
                                                      IO::LAS_property::X,
                                                      IO::LAS_property::Y,
                                                      IO::LAS_property::Z> point_property,
                                                      PropertyHandler&& ... properties)
{
  return IO::write_LAS_with_properties(os, points, point_property, std::forward<PropertyHandler>(properties)...);
}

/**
   \ingroup PkgPointSetProcessing3IODeprecated

  \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::write_LAS()` should be used instead.
*/
template <typename PointRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_las_points(std::ostream& os, const PointRange& points, const CGAL_BGL_NP_CLASS& np)
{
  return IO::write_LAS(os, points, np);
}

#endif //CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_POINT_SET_PROCESSING_WRITE_LAS_POINTS_H
