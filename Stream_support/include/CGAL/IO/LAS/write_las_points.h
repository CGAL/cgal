// Copyright (c) 2017 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_IO_LAS_WRITE_LAS_POINTS_H
#define CGAL_IO_LAS_WRITE_LAS_POINTS_H

#include <CGAL/IO/helpers.h>
#include <CGAL/IO/LAS/Las_property.h>
#include <CGAL/IO/LAS.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/assertions.h>

#include <boost/cstdint.hpp>
#include <boost/version.hpp>

#ifdef BOOST_MSVC
#  pragma warning(push)
#  pragma warning(disable:4251) // DLL warning from LASlib
#endif

#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

#define USE_AS_DLL 1
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
#include <type_traits>

namespace CGAL {

namespace IO {

// documented in ../LAS.h
template <typename PointMap>
std::tuple<PointMap, LAS_property::X, LAS_property::Y, LAS_property::Z >
make_las_point_writer(PointMap point_map)
{
  return std::make_tuple (point_map, LAS_property::X(), LAS_property::Y(), LAS_property::Z());
}

/// \cond SKIP_IN_MANUAL

namespace internal {
namespace LAS {

  inline void output_value(LASpoint& r, const unsigned short& v, const LAS_property::Intensity&)
  { r.set_intensity(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, const LAS_property::Return_number&)
  { r.set_return_number(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, const LAS_property::Number_of_returns&)
  { r.set_number_of_returns(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, const LAS_property::Scan_direction_flag&)
  { r.set_scan_direction_flag(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, const LAS_property::Edge_of_flight_line&)
  { r.set_edge_of_flight_line(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, const LAS_property::Classification&)
  { r.set_classification(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, const LAS_property::Synthetic_flag&)
  { r.set_synthetic_flag(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, const LAS_property::Keypoint_flag&)
  { r.set_keypoint_flag(v); }
  inline void output_value(LASpoint& r, const unsigned char& v, const LAS_property::Withheld_flag&)
  { r.set_withheld_flag(v); }
  inline void output_value(LASpoint& r, const float& v, const LAS_property::Scan_angle&)
  {
#if LAS_TOOLS_VERSION < 250517
    r.set_scan_angle_rank(I8_QUANTIZE(v));
#else
    r.set_scan_angle(v);
#endif
  }
  inline void output_value(LASpoint& r, const unsigned char& v, const LAS_property::User_data&)
  { r.set_user_data(v); }
  inline void output_value(LASpoint& r, const unsigned short& v, const LAS_property::Point_source_ID&)
  { r.set_point_source_ID(v); }
  inline void output_value(LASpoint& r, const unsigned int& v, const LAS_property::Deleted_flag&)
  { r.set_deleted_flag(v); }
  inline void output_value(LASpoint& r, const double& v, const LAS_property::GPS_time&)
  { r.set_gps_time(v); }
  inline void output_value(LASpoint& r, const unsigned short& v, const LAS_property::R&)
  { r.set_R(v); }
  inline void output_value(LASpoint& r, const unsigned short& v, const LAS_property::G&)
  { r.set_G(v); }
  inline void output_value(LASpoint& r, const unsigned short& v, const LAS_property::B&)
  { r.set_B(v); }
  inline void output_value(LASpoint& r, const unsigned short& v, const LAS_property::I&)
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

  template<typename Value, typename Tuple, std::size_t... Is>
  void output_tuple(LASpoint& point, const Value& v, const Tuple& t, std::index_sequence<Is...>) {
    (output_value(point, std::get<Is>(v), std::get<Is>(t)), ...);
  }

  template <typename ForwardIterator,
            typename PropertyMap,
            typename ... T>
  void output_properties(LASpoint& point,
                         ForwardIterator it,
                         std::tuple<PropertyMap, T ...>&& current)
  {
    output_tuple(point, get(std::get<0>(current), *it), std::tuple<T ...>(), std::index_sequence_for<T ...>{});
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

  template <typename ForwardIterator,
            typename PropertyMap,
            typename ... T,
            typename NextPropertyHandler,
            typename ... PropertyHandler>
  void output_properties(LASpoint& point,
                         ForwardIterator it,
                         std::tuple<PropertyMap, T ...>&& current,
                         NextPropertyHandler&& next,
                         PropertyHandler&& ... properties)
  {
    output_tuple(point, get(std::get<0>(current), *it), std::tuple<T ...>(), std::index_sequence_for<T ...>{});
    output_properties(point, it, std::forward<NextPropertyHandler>(next),
                      std::forward<PropertyHandler>(properties)...);
  }

} // namespace LAS
} // namespace internal

/// \endcond

// documented in ../LAS.h
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
  CGAL_precondition(points.begin() != points.end());

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

//  documented in ../LAS.h
template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
bool write_LAS(std::ostream& os,
               const PointRange& points,
               const CGAL_NP_CLASS& np,
               std::enable_if_t<internal::is_Range<PointRange>::value>*
               )
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetPointMap<PointRange, CGAL_NP_CLASS>::type PointMap;
  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));

  if(!os)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  return write_LAS_with_properties(os, points, make_las_point_writer(point_map));
}

// documented in ../LAS.h
template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
bool write_LAS(const std::string& filename,
               const PointRange& points,
               const CGAL_NP_CLASS& np,
               std::enable_if_t<internal::is_Range<PointRange>::value>*
               )
{
  std::ofstream os(filename, std::ios::binary);
  CGAL::IO::set_mode(os, CGAL::IO::BINARY);
  return write_LAS(os, points, np);
}

} // namespace IO

} // namespace CGAL

#endif // CGAL_IO_LAS_WRITE_LAS_POINTS_H
