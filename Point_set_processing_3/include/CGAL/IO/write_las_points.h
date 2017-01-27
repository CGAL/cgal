// Copyright (c) 2017  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s) : Simon Giraudot

#ifndef CGAL_WRITE_LAS_POINTS_H
#define CGAL_WRITE_LAS_POINTS_H

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Kernel_traits.h>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

#include <iostream>
#include <sstream>
#include <string>

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include <laswriter_las.hpp>
#pragma GCC diagnostic pop
#endif

namespace CGAL {


namespace LAS
{

  /**
     \ingroup PkgPointSetProcessing
     
     Generates a LAS property handler to write 3D points. 

     \sa `write_las_points_with_properties()`

     \tparam PointMap the property map used to store points.
  */
  template <typename PointMap>
  cpp11::tuple<PointMap, Property::X, Property::Y, Property::Z >
  make_point_writer(PointMap point_map)
  {
    return cpp11::make_tuple (point_map, Property::X(), Property::Y(), Property::Z());
  }

  /// \cond SKIP_IN_MANUAL
  
namespace internal {

  void output_value(LASpoint& r, unsigned short& v, LAS::Property::Intensity&)
  { r.set_intensity(v); }
  void output_value(LASpoint& r, unsigned char& v, LAS::Property::Return_number&)
  { r.set_return_number(v); }
  void output_value(LASpoint& r, unsigned char& v, LAS::Property::Number_of_returns&)
  { r.set_number_of_returns(v); }
  void output_value(LASpoint& r, unsigned char& v, LAS::Property::Scan_direction_flag&)
  { r.set_scan_direction_flag(v); }
  void output_value(LASpoint& r, unsigned char& v, LAS::Property::Edge_of_flight_line&)
  { r.set_edge_of_flight_line(v); }
  void output_value(LASpoint& r, unsigned char& v, LAS::Property::Classification&)
  { r.set_classification(v); }
  void output_value(LASpoint& r, unsigned char& v, LAS::Property::Synthetic_flag&)
  { r.set_synthetic_flag(v); }
  void output_value(LASpoint& r, unsigned char& v, LAS::Property::Keypoint_flag&)
  { r.set_keypoint_flag(v); }
  void output_value(LASpoint& r, unsigned char& v, LAS::Property::Withheld_flag&)
  { r.set_withheld_flag(v); }
  void output_value(LASpoint& r, float& v, LAS::Property::Scan_angle&)
  { r.set_scan_angle_rank(char(v)); }
  void output_value(LASpoint& r, unsigned char& v, LAS::Property::User_data&)
  { r.set_user_data(v); }
  void output_value(LASpoint& r, unsigned short& v, LAS::Property::Point_source_ID&)
  { r.set_point_source_ID(v); }
  void output_value(LASpoint& r, unsigned int& v, LAS::Property::Deleted_flag&)
  { r.set_deleted_flag(v); }
  void output_value(LASpoint& r, double& v, LAS::Property::GPS_time&)
  { r.set_gps_time(v); }
  void output_value(LASpoint& r, unsigned short& v, LAS::Property::R&)
  { r.set_R(v); }
  void output_value(LASpoint& r, unsigned short& v, LAS::Property::G&)
  { r.set_G(v); }
  void output_value(LASpoint& r, unsigned short& v, LAS::Property::B&)
  { r.set_B(v); }
  void output_value(LASpoint& r, unsigned short& v, LAS::Property::I&)
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
                          std::pair<PropertyMap, T>& current)
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
                          std::pair<PropertyMap, T>& current,
                          NextPropertyHandler& next,
                          PropertyHandler&& ... properties)
  {
    output_value (point, get(current.first, *it), current.second);
    output_properties (point, it, next, properties...);
  }

  
} // namespace internal
  

/// \endcond

}

//===================================================================================
/// \ingroup PkgPointSetProcessing
/// Saves the [first, beyond) range of points with properties to a
/// .las stream.
///
/// Properties are handled through a variadic list of property
/// handlers. A property handle is a `std::pair<PropertyMap,
/// LasProperty >` used to write a scalar value `LasProperty::type` as
/// a LAS property (for example, writing an `int` vairable as an `int`
/// LAS property). An exception is used for points that are written
/// using a `cpp11::tuple` object.
///
/// @sa `LAS::make_point_writer()`
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointMap is a model of `ReadablePropertyMap` with a value_type = `CGAL::Point_3`.
/// @tparam PropertyHandler handlers to recover properties.
///
/// @return true on success.
template <typename ForwardIterator,
          typename PointMap,
          typename ... PropertyHandler>
bool write_las_points_with_properties (std::ostream& stream,  ///< output stream.
                                       ForwardIterator first, ///< iterator over the first input point.
                                       ForwardIterator beyond, ///< past-the-end iterator over the input points.
                                       cpp11::tuple<PointMap,
                                                    LAS::Property::X,
                                                    LAS::Property::Y,
                                                    LAS::Property::Z> point_property,  ///< property handler for points
                                       PropertyHandler&& ... properties) ///< parameter pack of property handlers
{
  CGAL_point_set_processing_precondition(first != beyond);

  if(!stream)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  CGAL::Bbox_3 bbox = CGAL::bbox_3
    (boost::make_transform_iterator (first, CGAL::Property_map_to_unary_function<PointMap>(std::get<0>(point_property))),
     boost::make_transform_iterator (beyond, CGAL::Property_map_to_unary_function<PointMap>(std::get<0>(point_property))));
  
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
  for(ForwardIterator it = first; it != beyond; it++)
  {
    const typename PointMap::value_type& p = get(std::get<0>(point_property), *it);
    laspoint.set_X ((unsigned int)((p.x() - header.x_offset) / header.x_scale_factor));
    laspoint.set_Y ((unsigned int)((p.y() - header.y_offset) / header.y_scale_factor));
    laspoint.set_Z ((unsigned int)((p.z() - header.z_offset) / header.z_scale_factor));
    LAS::internal::output_properties (laspoint, it, properties...);

    laswriter.write_point (&laspoint);
    laswriter.update_inventory(&laspoint);
  }

  laswriter.update_header(&header, TRUE);

  laswriter.close();

  return ! stream.fail();
}

//===================================================================================
/// \ingroup PkgPointSetProcessing
/// Saves the [first, beyond) range of points (positions only) to a
/// .las stream. 
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointMap is a model of `ReadablePropertyMap` with a value_type = `Point_3<Kernel>`.
///        It can be omitted if the value type of `ForwardIterator` is convertible to `Point_3<Kernel>`.
///
/// @return true on success.

// This variant requires all parameters.
template < typename ForwardIterator,
           typename PointMap >
bool
write_las_points(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond, ///< past-the-end input point.
  PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return write_las_points_with_properties (stream, first, beyond, LAS::make_point_writer(point_map));
}

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
template < typename ForwardIterator >
bool
write_las_points(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond) ///< past-the-end input point.
{
  return write_las_points(
    stream,
    first, beyond,
    make_identity_property_map(
    typename std::iterator_traits<ForwardIterator>::value_type())
    );
}
/// @endcond


} //namespace CGAL

#endif // CGAL_WRITE_LAS_POINTS_H
