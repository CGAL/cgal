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

#ifndef CGAL_READ_PLY_POINTS_H
#define CGAL_READ_PLY_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/config.h>
#if defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE) || defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES)
#error CGAL PLY reader requires a C++11 compiler
#endif

#include <tuple>

#include <CGAL/IO/PLY.h>
#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/IO/io.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

#include <iostream>
#include <sstream>
#include <string>

namespace CGAL {

#ifdef DOXYGEN_RUNNING // Document some parts from Stream_support here for convenience
  /**
     \ingroup PkgPointSetProcessing3IOPly
     
     Class used to identify a %PLY property as a type and a name.

     \sa `read_ply_points_with_properties()`
  */
  template <typename T>
  struct PLY_property
  {
    typedef T type;
    const char* name;
    PLY_property (const char* name) : name (name) { }
  };

  /**
     \ingroup PkgPointSetProcessing3IOPly
     
     Generates a %PLY property handler to read 3D points. Points are
     constructed from the input using 3 %PLY properties of type `FT`
     and named `x`, `y` and `z`. `FT` is `float` if the points use
     `CGAL::Simple_cartesian<float>` and `double` otherwise.

     \sa `read_ply_points_with_properties()`

     \tparam PointMap the property map used to store points.
  */
  template <typename PointMap>
  std::tuple<PointMap,
             typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3,
             PLY_property<FT>, PLY_property<FT>, PLY_property<FT> >
  make_ply_point_reader(PointMap point_map);

  /**
     \ingroup PkgPointSetProcessing3IOPly
     
     Generates a %PLY property handler to read 3D normal
     vectors. Vectors are constructed from the input using 3 PLY
     properties of type `FT` and named `nx`, `ny` and `nz`. `FT`
     is `float` if the points use `CGAL::Simple_cartesian<float>` and
     `double` otherwise.

     \sa `read_ply_points_with_properties()`

     \tparam VectorMap the property map used to store vectors.
  */
  template <typename VectorMap>
  std::tuple<VectorMap,
             typename Kernel_traits<typename VectorMap::value_type>::Kernel::Construct_vector_3,
             PLY_property<FT>, PLY_property<FT>, PLY_property<FT> >
  make_ply_normal_reader(VectorMap normal_map);
#endif // DOXYGEN_RUNNING

/*
  \ingroup PkgPointSetProcessing3IOPly

  Reads user-selected points properties from a .ply stream (ASCII or
  binary).
  Potential additional point properties and faces are ignored.

  Properties are handled through a variadic list of property
  handlers. A `PropertyHandler` can either be:

  - A `std::pair<PropertyMap, PLY_property<T> >` if the user wants
  to read a %PLY property as a scalar value T (for example, storing
  an `int` %PLY property into an `int` variable).

  - A `std::tuple<PropertyMap, Constructor,
  PLY_property<T>...>` if the user wants to use one or several PLY
  properties to construct a complex object (for example, storing 3
  `uchar` %PLY properties into a %Color object that can for example
  be a `CGAL::cpp11::array<unsigned char, 3>`). In that case, the
  second element of the tuple should be a functor that constructs
  the value type of `PropertyMap` from N objects of types `T`.

  \sa `make_ply_point_reader()`
  \sa `make_ply_normal_reader()`

  \cgalRequiresCPP11

  \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
  It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
  \tparam OutputIterator iterator over output points.
  \tparam PropertyHandler handlers to recover properties.

  \return `true` on success.
*/
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename ... PropertyHandler>
bool read_ply_points_with_properties (std::istream& stream,
                                      OutputIterator output,
                                      PropertyHandler&& ... properties)
{
  typedef typename value_type_traits<OutputIterator>::type OutputValueType;
    
  if(!stream)
    {
      std::cerr << "Error: cannot open file" << std::endl;
      return false;
    }

  internal::PLY::PLY_reader reader;
  
  if (!(reader.init (stream)))
  {
    stream.setstate(std::ios::failbit);
    return false;
  }
  
  for (std::size_t i = 0; i < reader.number_of_elements(); ++ i)
  {
    internal::PLY::PLY_element& element = reader.element(i);

    for (std::size_t j = 0; j < element.number_of_items(); ++ j)
    {
      for (std::size_t k = 0; k < element.number_of_properties(); ++ k)
      {
        internal::PLY::PLY_read_number* property = element.property(k);
        property->get (stream);

        if (stream.fail())
          return false;
      }

      if (element.name() == "vertex" || element.name() == "vertices")
      {
        OutputValueType new_element;
        internal::PLY::process_properties (element, new_element, std::forward<PropertyHandler>(properties)...);
        *(output ++) = new_element;
      }
    }
  }

  return true;
}

/// \cond SKIP_IN_MANUAL
template <typename OutputIterator,
          typename ... PropertyHandler>
bool read_ply_points_with_properties (std::istream& stream,
                                      OutputIterator output,
                                      PropertyHandler&& ... properties)
{
  typedef typename value_type_traits<OutputIterator>::type OutputValueType;
 
  return read_ply_points_with_properties<OutputValueType>
    (stream, output, std::forward<PropertyHandler>(properties)...);
}
/// \endcond
  
/**
   \ingroup PkgPointSetProcessing3IOPly
   Reads points (positions + normals, if available) from a .ply
   stream (ASCII or binary).
   Potential additional point properties and faces are ignored.

   \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
   It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
   \tparam OutputIterator iterator over output points.

   \param stream input stream.
   \param output output iterator over points.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `WritablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `ReadWritePropertyMap` with value type
     `geom_traits::Vector_3`. If this parameter is omitted, normals in the input stream are
     ignored.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

   \return true on success.

   \cgalRequiresCPP11
*/
template < typename OutputIteratorValueType,
           typename OutputIterator,
#ifdef DOXYGEN_RUNNING
           typename NamedParameters
#else
           typename CGAL_BGL_NP_TEMPLATE_PARAMETERS
#endif
>
bool read_ply_points(std::istream& stream,
                     OutputIterator output,
#ifdef DOXYGEN_RUNNING
                     const NamedParameters& np)
#else
                     const CGAL_BGL_NP_CLASS& np)
#endif
{
  using boost::choose_param;

  typedef Point_set_processing_3::Fake_point_range<OutputIteratorValueType> PointRange;
  
  // basic geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, CGAL_BGL_NP_CLASS>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, CGAL_BGL_NP_CLASS>::type NormalMap;

  bool has_normals = !(boost::is_same<NormalMap,
                       typename Point_set_processing_3::GetNormalMap<PointRange, CGAL_BGL_NP_CLASS>::NoMap>::value);

  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
  NormalMap normal_map = choose_param(get_param(np, internal_np::normal_map), NormalMap());

  if (has_normals)
    return read_ply_points_with_properties (stream, output,
                                            make_ply_point_reader (point_map),
                                            make_ply_normal_reader (normal_map));
  // else
  return read_ply_points_with_properties (stream, output,
                                          make_ply_point_reader (point_map));
}


/// \cond SKIP_IN_MANUAL  
// variant with default NP
template <typename OutputIteratorValueType,
          typename OutputIterator>
bool
read_ply_points(
  std::istream& stream, ///< input stream.
  OutputIterator output) ///< output iterator over points.
{
  return read_ply_points<OutputIteratorValueType>
    (stream, output, CGAL::parameters::all_default());
}

// variant with default output iterator value type
template <typename OutputIterator,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool
read_ply_points(
  std::istream& stream, ///< input stream.
  OutputIterator output,
  const CGAL_BGL_NP_CLASS& np)
{
  return read_ply_points<typename value_type_traits<OutputIterator>::type>
    (stream, output, np);
}

// variant with default NP and output iterator value type
template <typename OutputIterator>
bool
read_ply_points(
  std::istream& stream, ///< input stream.
  OutputIterator output)
{
  return read_ply_points<typename value_type_traits<OutputIterator>::type>
    (stream, output, CGAL::parameters::all_default());
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointMap,
           typename NormalMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_ply_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map));
}

// deprecated API
template < typename OutputIterator,
           typename PointMap,
           typename NormalMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_ply_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map));
}

// deprecated API
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename NormalMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_ply_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::normal_map (normal_map));
}

// deprecated API
template < typename OutputIterator,
           typename NormalMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_ply_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::normal_map (normal_map));
}

// deprecated API  
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points(), please update your code")
bool
read_ply_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_ply_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::point_map (point_map));
}

// deprecated API
template < typename OutputIterator,
           typename PointMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points(), please update your code")
bool read_ply_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_ply_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::point_map (point_map));
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond


} //namespace CGAL

#undef TRY_TO_GENERATE_POINT_PROPERTY
#undef TRY_TO_GENERATE_SIZED_FACE_PROPERTY
#undef TRY_TO_GENERATE_FACE_PROPERTY

#endif // CGAL_READ_PLY_POINTS_H
