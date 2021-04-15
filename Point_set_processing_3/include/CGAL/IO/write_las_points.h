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

#ifndef CGAL_WRITE_LAS_POINTS_H
#define CGAL_WRITE_LAS_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/config.h>

#include <tuple>

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Kernel_traits.h>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

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
#include <sstream>
#include <string>


namespace CGAL {


  /**
     \ingroup PkgPointSetProcessing3IOLas

     Generates a %LAS property handler to write 3D points.

     \sa `write_las_points_with_properties()`

     \tparam PointMap the property map used to store points.
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
  void output_properties (LASpoint&,
                          ForwardIterator)
  {
  }

  template <typename ForwardIterator,
            typename PropertyMap,
            typename T>
  void output_properties (LASpoint& point,
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
  void output_properties (LASpoint& point,
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
   Saves the range of `points` with properties to a
   .las stream.

   Properties are handled through a variadic list of property
   handlers. A `PropertyHandle` is a `std::pair<PropertyMap,
   LAS_property::Tag >` used to write a scalar value
   `LAS_property::Tag::type` as a %LAS property (for example,
   writing an `int` vairable as an `int` %LAS property). An exception
   is used for points that are written using a `std::tuple` object.

   See documentation of `read_las_points_with_properties()` for the
   list of available `LAS_property::Tag` classes.

   \sa `make_las_point_writer()`

   \cgalRequiresCPP11

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam PointMap is a model of `ReadablePropertyMap` with a value_type = `CGAL::Point_3`.
   \tparam PropertyHandler handlers to recover properties.

   \return `true` on success.
*/
template <typename PointRange,
          typename PointMap,
          typename ... PropertyHandler>
bool write_las_points_with_properties (std::ostream& stream,  ///< output stream.
                                       const PointRange& points, ///< input point range.
                                       std::tuple<PointMap,
                                       LAS_property::X,
                                       LAS_property::Y,
                                       LAS_property::Z> point_property,  ///< property handler for points
                                       PropertyHandler&& ... properties) ///< parameter pack of property handlers
{
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  if(!stream)
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
  laswriter.open (stream, &header);

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

  return ! stream.fail();
}

/**
   \ingroup PkgPointSetProcessing3IOLas
   Saves the range of `points` (positions only) to a
   .las stream.

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
     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \return true on success.
   \cgalRequiresCPP11
*/
template < typename PointRange,
           typename NamedParameters>
bool
write_las_points(
  std::ostream& stream,
  const PointRange& points,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::type PointMap;
  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));

  return write_las_points_with_properties (stream, points, make_las_point_writer(point_map));
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template < typename PointRange>
bool
write_las_points(
  std::ostream& stream,
  const PointRange& points)
{
  return write_las_points
    (stream, points, CGAL::Point_set_processing_3::parameters::all_default(points));
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template < typename ForwardIterator,
           typename PointMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_las_points(), please update your code")
bool
write_las_points(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond, ///< past-the-end input point.
  PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_las_points
    (stream, points,
     CGAL::parameters::point_map(point_map));
}

// deprecated API
template < typename ForwardIterator >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_las_points(), please update your code")
bool
write_las_points(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond) ///< past-the-end input point.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_las_points
    (stream, points);
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond


} //namespace CGAL

#endif // CGAL_WRITE_LAS_POINTS_H
