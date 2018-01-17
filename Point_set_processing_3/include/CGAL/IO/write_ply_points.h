// Copyright (c) 2015  Geometry Factory
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Simon Giraudot

#ifndef CGAL_WRITE_PLY_POINTS_H
#define CGAL_WRITE_PLY_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/config.h>
#if defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE) || defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES)
#error CGAL PLY writer requires a C++11 compiler
#endif

#include <tuple>

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/Iterator_range.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/version.hpp>

#include <iostream>
#include <iterator>

namespace CGAL {

  /**
     \ingroup PkgPointSetProcessingIOPly
     
     Generates a %PLY property handler to write 3D points. Points are
     written as 3 %PLY properties of type `double` and named `x`, `y`
     and `z`.

     \sa `write_ply_points_with_properties()`

     \tparam PointMap the property map used to store points.
  */
  template <typename PointMap>
  std::tuple<PointMap, PLY_property<double>, PLY_property<double>, PLY_property<double> >
  make_ply_point_writer(PointMap point_map)
  {
    return std::make_tuple (point_map, PLY_property<double>("x"), PLY_property<double>("y"), PLY_property<double>("z"));
  }

  /**
     \ingroup PkgPointSetProcessingIOPly
     
     Generates a %PLY property handler to write 3D normal
     vectors. Vectors are written as 3 %PLY properties of type `double`
     and named `nx`, `ny` and `nz`.

     \sa `write_ply_points_with_properties()`

     \tparam VectorMap the property map used to store vectors.
  */
  template <typename VectorMap>
  std::tuple<VectorMap, PLY_property<double>, PLY_property<double>, PLY_property<double> >
  make_ply_normal_writer(VectorMap normal_map)
  {
    return std::make_tuple (normal_map, PLY_property<double>("nx"), PLY_property<double>("ny"), PLY_property<double>("nz"));
  }

  /// \cond SKIP_IN_MANUAL

namespace internal {

  namespace PLY {

  template <typename T> void property_header_type (std::ostream& stream) { stream << "undefined_type"; }

  template <> void property_header_type<char> (std::ostream& stream) { stream << "char"; }
  template <> void property_header_type<unsigned char> (std::ostream& stream) { stream << "uchar"; }
  template <> void property_header_type<short> (std::ostream& stream) { stream << "short"; }
  template <> void property_header_type<unsigned short> (std::ostream& stream) { stream << "ushort"; }
  template <> void property_header_type<int> (std::ostream& stream) { stream << "int"; }
  template <> void property_header_type<unsigned int> (std::ostream& stream) { stream << "uint"; }
  template <> void property_header_type<float> (std::ostream& stream) { stream << "float"; }
  template <> void property_header_type<double> (std::ostream& stream) { stream << "double"; }
    
  template <typename T>
  void property_header (std::ostream& stream, const PLY_property<T>& prop)
  {
    stream << "property ";
    property_header_type<T>(stream);
    stream << " " << prop.name << std::endl;
  }

  template <typename T>
  void property_header (std::ostream& stream, const PLY_property<std::vector<T> >& prop)
  {
    stream << "property list uchar ";
    property_header_type<T>(stream);
    stream << " " << prop.name << std::endl;
  }

  
  template <std::size_t N>
  struct Properties_header
  {
    template <class PLY_property_tuple>
    static void write(std::ostream& stream, PLY_property_tuple& wrappers)
    {
      Properties_header<N-1>::write(stream, wrappers);
      property_header (stream, std::get<N+1>(wrappers));
    }
  };
  template <>
  struct Properties_header<0>
  {
    template <class PLY_property_tuple>
    static void write(std::ostream& stream, PLY_property_tuple& wrappers)
    {
      property_header (stream, std::get<1>(wrappers));
    }
  };

  template <typename PropertyMap,
            typename ... T>
  void output_property_header (std::ostream& stream,
                               std::tuple<PropertyMap, PLY_property<T>... >&& current)
  {
    Properties_header<sizeof...(T)-1>::write(stream, current); 
  }


  template <typename PropertyMap,
            typename T>
  void output_property_header (std::ostream& stream,
                               std::pair<PropertyMap, PLY_property<T> >&& current)
  {
    property_header (stream, current.second);
  }

  template <typename PropertyMap,
            typename T,
            typename NextPropertyHandler,
            typename ... PropertyHandler>
  void output_property_header (std::ostream& stream,
                               std::pair<PropertyMap, PLY_property<T> >&& current,
                               NextPropertyHandler&& next,
                               PropertyHandler&& ... properties)
  {
    property_header (stream, current.second);
    output_property_header (stream, std::forward<NextPropertyHandler>(next),
                            std::forward<PropertyHandler>(properties)...);
  }
  template <typename PropertyMap,
            typename ... T,
            typename NextPropertyHandler,
            typename ... PropertyHandler>
  void output_property_header (std::ostream& stream,
                               std::tuple<PropertyMap, PLY_property<T>... >&& current,
                               NextPropertyHandler&& next,
                               PropertyHandler&& ... properties)
  {
    Properties_header<sizeof...(T)-1>::write(stream, current); 
    output_property_header (stream, std::forward<NextPropertyHandler>(next),
                            std::forward<PropertyHandler>(properties)...);
  }


  template <typename ForwardIterator,
            typename PropertyMap>
  void property_write (std::ostream& stream, ForwardIterator it, PropertyMap map)
  {
    stream << CGAL::oformat(get (map, *it));
  }

  template <typename ForwardIterator,
            typename PropertyMap,
            typename T>
  void simple_property_write (std::ostream& stream, ForwardIterator it,
                              std::pair<PropertyMap, PLY_property<T> >&& map)
  {
    if (CGAL::get_mode(stream) == IO::ASCII)
      stream << get (map.first, *it);
    else
      {
        typename PropertyMap::value_type value = get(map.first, *it);
        stream.write (reinterpret_cast<char*>(&value), sizeof(value));
      }
  }

  template <typename ForwardIterator,
            typename PropertyMap,
            typename T>
  void simple_property_write (std::ostream& stream, ForwardIterator it,
                              std::pair<PropertyMap, PLY_property<std::vector<T> > > map)
  {
    const typename PropertyMap::reference value = get(map.first, *it);
    
    if (CGAL::get_mode(stream) == IO::ASCII)
    {
      stream << value.size();
      for (std::size_t i = 0; i < value.size(); ++ i)
        stream << " " << value[i];
    }
    else
      {
        unsigned char size = static_cast<unsigned char>(value.size());
        stream.write (reinterpret_cast<char*>(&size), sizeof(size));
        for (std::size_t i = 0; i < value.size(); ++ i)
        {
          T t = T(value[i]);
          stream.write (reinterpret_cast<char*>(&t), sizeof(t));
        }
      }
  }

    
  template <typename ForwardIterator,
            typename PropertyMap,
            typename ... T>
  void output_properties (std::ostream& stream,
                          ForwardIterator it,
                          std::tuple<PropertyMap, PLY_property<T>... >&& current)
  {
    property_write (stream, it, std::get<0>(current));
    if (get_mode(stream) == IO::ASCII)
      stream << std::endl;
  }


  template <typename ForwardIterator,
            typename PropertyMap,
            typename T>
  void output_properties (std::ostream& stream,
                          ForwardIterator it,
                          std::pair<PropertyMap, PLY_property<T> >&& current)
  {
    simple_property_write (stream, it, std::forward<std::pair<PropertyMap, PLY_property<T> > >(current));
    if (get_mode(stream) == IO::ASCII)
      stream << std::endl;
  }

  template <typename ForwardIterator,
            typename PropertyMap,
            typename T,
            typename NextPropertyHandler,
            typename ... PropertyHandler>
  void output_properties (std::ostream& stream,
                          ForwardIterator it,
                          std::pair<PropertyMap, PLY_property<T> >&& current,
                          NextPropertyHandler& next,
                          PropertyHandler&& ... properties)
  {
    simple_property_write (stream, it, current);
    if (get_mode(stream) == IO::ASCII)
      stream << " ";
    output_properties (stream, it, std::forward<NextPropertyHandler>(next),
                       std::forward<PropertyHandler>(properties)...);
  }
  
  template <typename ForwardIterator,
            typename PropertyMap,
            typename ... T,
            typename NextPropertyHandler,
            typename ... PropertyHandler>
  void output_properties (std::ostream& stream,
                          ForwardIterator it,
                          std::tuple<PropertyMap, PLY_property<T>... >&& current,
                          NextPropertyHandler&& next,
                          PropertyHandler&& ... properties)
  {
    property_write (stream, it, std::get<0>(current));
    if (get_mode(stream) == IO::ASCII)
      stream << " ";
    output_properties (stream, it, std::forward<NextPropertyHandler>(next),
                       std::forward<PropertyHandler>(properties)...);
  }

  } // namespace PLY
    
} // namespace internal

  /// \endcond


/**
   \ingroup PkgPointSetProcessingIOPly
   Saves the range of `points` with properties to a
   .ply stream. %PLY is either ASCII or binary depending on the value
   of `CGAL::get_mode(stream)`.

   Properties are handled through a variadic list of property
   handlers. A `PropertyHandler` can either be:

   - A `std::pair<PropertyMap, PLY_property<T> >` if the user wants
   to write a scalar value T as a %PLY property (for example, writing
   an `int` variable as an `int` %PLY property).

   - A `std::tuple<PropertyMap, PLY_property<T>...>` if the
   user wants to write a complex object as several %PLY
   properties. In that case, a specialization of `Output_rep` must
   be provided for `PropertyMap::value_type` that handles both ASCII
   and binary output (see `CGAL::get_mode()`).

   \sa `make_ply_point_writer()`
   \sa `make_ply_normal_writer()`

   \cgalRequiresCPP11

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam PropertyHandler handlers to recover properties.

   \return `true` on success.
*/
template < typename PointRange,
           typename ... PropertyHandler>
bool
write_ply_points_with_properties(
  std::ostream& stream, ///< output stream.
  const PointRange& points, ///< input point range.
  PropertyHandler&& ... properties) ///< parameter pack of property handlers
{
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  if(!stream)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // Write header
  stream << "ply" << std::endl
         << ((get_mode(stream) == IO::BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
         << "comment Generated by the CGAL library" << std::endl
         << "element vertex " << points.size() << std::endl;
  
  internal::PLY::output_property_header (stream, std::forward<PropertyHandler>(properties)...);
  
  stream << "end_header" << std::endl;
  

  // Write positions + normals
  for(typename PointRange::const_iterator it = points.begin(); it != points.end(); it++)
  {
    internal::PLY::output_properties (stream, it, std::forward<PropertyHandler>(properties)...);
  }

  return ! stream.fail();
}

/**
   \ingroup PkgPointSetProcessingIOPly
   Saves the range of `points` (positions + normals, if available) to
   a .ply stream. %PLY is either ASCII or binary depending on the
   value of `CGAL::get_mode(stream)`.

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param stream output stream.
   \param points input point range.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `ReadablePropertyMap` with value type
     `geom_traits::Vector_3`.\cgalParamEnd If this parameter is omitted, normals are not written to the
     output stream.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

   \return true on success.
   \cgalRequiresCPP11
*/
template <typename PointRange,
          typename NamedParameters>
bool
write_ply_points(
  std::ostream& stream,
  const PointRange& points,
  const NamedParameters& np)
{
  using boost::choose_param;

  // basic geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;

  bool has_normals = !(boost::is_same<NormalMap,
                       typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value);

  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
  NormalMap normal_map = choose_param(get_param(np, internal_np::normal_map), NormalMap());
  
  if (has_normals)
    return write_ply_points_with_properties(
      stream, points,
      make_ply_point_writer(point_map),
      make_ply_normal_writer(normal_map));
  // else
  return write_ply_points_with_properties(
    stream, points,
    make_ply_point_writer(point_map));
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename PointRange>
bool
write_ply_points(
  std::ostream& stream,
  const PointRange& points)
{
  return write_ply_points
    (stream, points, CGAL::Point_set_processing_3::parameters::all_default(points));
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template < typename ForwardIterator,
           typename PointMap,
           typename VectorMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_ply_points_and_normals(), please update your code")
bool
write_ply_points_and_normals(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond, ///< past-the-end input point.
  PointMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
  VectorMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_ply_points
    (stream, points,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map));
}
  
// deprecated API
template <typename ForwardIterator,
          typename VectorMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_ply_points_and_normals(), please update your code")
bool
write_ply_points_and_normals(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond, ///< past-the-end input point.
  VectorMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_ply_points
    (stream, points,
     CGAL::parameters::normal_map (normal_map));
}

// deprecated API
template < typename ForwardIterator,
           typename PointMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_ply_points(), please update your code")
bool
write_ply_points(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond, ///< past-the-end input point.
  PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_ply_points
    (stream, points,
     CGAL::parameters::point_map(point_map));
}

// deprecated API
template < typename ForwardIterator >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_ply_points(), please update your code")
bool
write_ply_points(
  std::ostream& stream, ///< output stream.
  ForwardIterator first, ///< first input point.
  ForwardIterator beyond) ///< past-the-end input point.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return write_ply_points
    (stream, points);
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond


} //namespace CGAL

#endif // CGAL_WRITE_PLY_POINTS_H
