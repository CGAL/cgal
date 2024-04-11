// Copyright (c) 2015  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_POINT_SET_PROCESSING_READ_PLY_POINTS_H
#define CGAL_POINT_SET_PROCESSING_READ_PLY_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/config.h>

#include <CGAL/IO/PLY.h>
#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/IO/io.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>

namespace CGAL {

namespace IO {

#ifdef DOXYGEN_RUNNING // Document some parts from Stream_support here for convenience
/**
   \ingroup PkgPointSetProcessing3IOPly

   Class used to identify a %PLY property as a type and a name.

   \sa `read_PLY_with_properties()`
*/
template <typename T>
struct PLY_property
{
  typedef T type;
  const char* name;
  PLY_property(const char* name) : name(name) { }
};

/**
   \ingroup PkgPointSetProcessing3IOPly

   Generates a %PLY property handler to read 3D points. Points are
   constructed from the input using 3 %PLY properties of type `FT`
   and named `x`, `y` and `z`. `FT` is `float` if the points use
   `CGAL::Simple_cartesian<float>` and `double` otherwise.

   \tparam PointMap the property map used to store points.

   \sa `read_PLY_with_properties()`
   \sa \ref IOStreamPLY
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

   \tparam VectorMap the property map used to store vectors.

   \sa `read_PLY_with_properties()`
   \sa \ref IOStreamPLY
*/
template <typename VectorMap>
std::tuple<VectorMap,
           typename Kernel_traits<typename VectorMap::value_type>::Kernel::Construct_vector_3,
           PLY_property<FT>, PLY_property<FT>, PLY_property<FT> >
make_ply_normal_reader(VectorMap normal_map);

#endif // DOXYGEN_RUNNING

/**
  \ingroup PkgPointSetProcessing3IOPly

  \brief reads user-selected points properties from a .ply stream (ASCII or binary).

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
  be a `std::array<unsigned char, 3>`). In that case, the
  second element of the tuple should be a functor that constructs
  the value type of `PropertyMap` from N objects of types `T`.

  \attention To read a binary file, the flag `std::ios::binary` must be set during the creation of the `ifstream`.

  \tparam OutputIteratorValueType type of objects that can be put in `PointOutputIterator`.
  It must be a model of `DefaultConstructible` and defaults to `value_type_traits<PointOutputIterator>::%type`.
  It can be omitted when the default is fine.
  \tparam PointOutputIterator iterator over output points.
  \tparam PropertyHandler handlers to recover properties.

  \returns `true` if reading was successful, `false` otherwise.

  \sa \ref IOStreamPLY
  \sa `make_ply_point_reader()`
  \sa `make_ply_normal_reader()`
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename ... PropertyHandler>
bool read_PLY_with_properties(std::istream& is,
                              PointOutputIterator output,
                              PropertyHandler&& ... properties)
{
  typedef typename value_type_traits<PointOutputIterator>::type OutputValueType;

  if(!is)
    return false;

  internal::PLY_reader reader(true);

  if(!(reader.init(is)))
  {
    is.setstate(std::ios::failbit);
    return false;
  }

  for(std::size_t i = 0; i < reader.number_of_elements(); ++ i)
  {
    internal::PLY_element& element = reader.element(i);

    for(std::size_t j = 0; j < element.number_of_items(); ++ j)
    {
      for(std::size_t k = 0; k < element.number_of_properties(); ++ k)
      {
        internal::PLY_read_number* property = element.property(k);
        property->get(is);

        if(is.fail())
          return false;
      }

      if(element.name() == "vertex" || element.name() == "vertices")
      {
        OutputValueType new_element;
        internal::process_properties(element, new_element, std::forward<PropertyHandler>(properties)...);
        *(output ++) = new_element;
      }
    }
  }

  return true;
}

/// \cond SKIP_IN_MANUAL

template <typename OutputIterator,
          typename ... PropertyHandler>
bool read_PLY_with_properties(std::istream& is,
                              OutputIterator output,
                              PropertyHandler&& ... properties)
{
  typedef typename value_type_traits<OutputIterator>::type OutputValueType;

  return read_PLY_with_properties<OutputValueType>(is, output, std::forward<PropertyHandler>(properties)...);
}

/// \endcond

/**
   \ingroup PkgPointSetProcessing3IOPly

   \brief reads points (positions + normals, if available), using the \ref IOStreamPLY.

   Potential additional point properties and faces are ignored.

  \attention To read a binary file, the flag `std::ios::binary` must be set during the creation of the `ifstream`.

   \tparam OutputIteratorValueType type of objects that can be put in `PointOutputIterator`.
   It must be a model of `DefaultConstructible` and defaults to `value_type_traits<PointOutputIterator>::%type`.
   It can be omitted when the default is fine.
   \tparam PointOutputIterator iterator over output points.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param is input stream.
   \param output output iterator over points.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals in the input stream are ignored.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if reading was successful, `false` otherwise.

   \sa `read_PLY_with_properties()`
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(std::istream& is,
              PointOutputIterator output,
              const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
              , std::enable_if_t<CGAL::is_iterator<PointOutputIterator>::value>* = nullptr
#endif
              )
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3::Fake_point_range<OutputIteratorValueType> PointRange;

  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, CGAL_NP_CLASS> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;

  PointMap point_map = NP_helper::get_point_map(np);
  NormalMap normal_map = NP_helper::get_normal_map(np);

  return read_PLY_with_properties(is, output,
                                  make_ply_point_reader(point_map),
                                  make_ply_normal_reader(normal_map));
}

/**
   \ingroup PkgPointSetProcessing3IOPly

   \brief reads points (positions + normals, if available), using the \ref IOStreamPLY.

   Potential additional point properties and faces are ignored.

   \tparam OutputIteratorValueType type of objects that can be put in `PointOutputIterator`.
   It must be a model of `DefaultConstructible` and defaults to `value_type_traits<PointOutputIterator>::%type`.
   It can be omitted when the default is fine.
   \tparam PointOutputIterator iterator over output points.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param fname input file name.
   \param output output iterator over points.
   \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamNBegin{use_binary_mode}
       \cgalParamDescription{indicates whether data should be read in binary (`true`) or in \ascii (`false`)}
       \cgalParamType{Boolean}
       \cgalParamDefault{`true`}
     \cgalParamNEnd

     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals in the input stream are ignored.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if reading was successful, `false` otherwise.

   \sa \ref IOStreamPLY
   \sa `read_PLY_with_properties()`
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(const std::string& fname,
              PointOutputIterator output,
              const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
              , std::enable_if_t<CGAL::is_iterator<PointOutputIterator>::value>* = nullptr
#endif
              )
{
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ifstream is(fname, std::ios::binary);
    CGAL::IO::set_mode(is, CGAL::IO::BINARY);
    return read_PLY<OutputIteratorValueType>(is, output, np);
  }
  else
  {
    std::ifstream is(fname);
    CGAL::IO::set_mode(is, CGAL::IO::ASCII);
    return read_PLY<OutputIteratorValueType>(is, output, np);
  }
}

/// \cond SKIP_IN_MANUAL

// variants with default output iterator value type
template <typename OutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(std::istream& is, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values(),
              std::enable_if_t<CGAL::is_iterator<OutputIterator>::value>* = nullptr)
{
  return read_PLY<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

template <typename OutputIterator,typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(const std::string& fname, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values(),
              std::enable_if_t<CGAL::is_iterator<OutputIterator>::value>* = nullptr)
{
  return read_PLY<typename value_type_traits<OutputIterator>::type>(fname, output, np);
}

/// \endcond

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/// \cond SKIP_IN_MANUAL

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointMap,
          typename NormalMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_PLY<OutputIteratorValueType>(is, output, parameters::point_map(point_map)
                                                                      .normal_map(normal_map));
}

template <typename OutputIterator,
          typename PointMap,
          typename NormalMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_PLY<typename value_type_traits<OutputIterator>::type>(is, output, parameters::point_map(point_map)
                                                                                               .normal_map(normal_map));
}

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename NormalMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_PLY<OutputIteratorValueType>(is, output, parameters::normal_map(normal_map));
}

template <typename OutputIterator, typename NormalMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points_and_normals(), please update your code")
bool read_ply_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_PLY<typename value_type_traits<OutputIterator>::type>(is, output, parameters::normal_map(normal_map));
}

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points(), please update your code")
bool read_ply_points(std::istream& is, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return IO::read_PLY<OutputIteratorValueType>(is, output, parameters::point_map(point_map));
}

template <typename OutputIterator,
          typename PointMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_ply_points(), please update your code")
bool read_ply_points(std::istream& is, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return IO::read_PLY<typename value_type_traits<OutputIterator>::type>(is, output, parameters::point_map(point_map));
}

/// \endcond

template <typename OutputIteratorValueType, typename OutputIterator, typename ... PropertyHandler>
CGAL_DEPRECATED bool read_ply_points_with_properties(std::istream& is, OutputIterator output, PropertyHandler&& ... properties)
{
  return IO::read_PLY_with_properties(is, output, std::forward<PropertyHandler>(properties)...);
}


template <typename OutputIteratorValueType, typename OutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool read_ply_points(std::istream& is, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values())
{
  return IO::read_PLY(is, output, np);
}

/// \cond SKIP_IN_MANUAL

template <typename OutputIterator,
          typename ... PropertyHandler>
CGAL_DEPRECATED bool read_ply_points_with_properties(std::istream& is, OutputIterator output, PropertyHandler&& ... properties)
{
  return IO::read_PLY_with_properties<typename value_type_traits<OutputIterator>::type>(is, output, std::forward<PropertyHandler>(properties)...);
}

// variant with default output iterator value type
template <typename OutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool read_ply_points(std::istream& is, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values())
{
  return IO::read_PLY<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

/// \endcond

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#undef TRY_TO_GENERATE_POINT_PROPERTY
#undef TRY_TO_GENERATE_SIZED_FACE_PROPERTY
#undef TRY_TO_GENERATE_FACE_PROPERTY

#endif // CGAL_POINT_SET_PROCESSING_READ_PLY_POINTS_H
