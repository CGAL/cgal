// Copyright (c) 2017 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_IO_LAS_H
#define CGAL_IO_LAS_H

#include <fstream>
#include <iostream>
#include <vector>
#include <type_traits>
#include <CGAL/IO/LAS/Las_property.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {

/**
   \ingroup PkgStreamSupportIoFuncsLAS

   generates a %LAS property handler to read 3D points. Points are
   constructed from the input the using 3 %LAS properties
   `LAS_property::X`, `LAS_property::Y` and `LAS_property::Z`.

   \tparam PointMap the property map used to store points.

   \sa `read_LAS_with_properties()`
   \sa \ref IOStreamLAS
*/
template <typename PointMap>
std::tuple<PointMap,
           typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3,
           LAS_property::X, LAS_property::Y, LAS_property::Z >
make_las_point_reader(PointMap point_map);


/**
   \ingroup PkgStreamSupportIoFuncsLAS

   \brief reads user-selected points properties from a .las or .laz stream.
   Potential additional properties are ignored.

   Properties are handled through a variadic list of property
   handlers. A `PropertyHandler` can either be:

   - A `std::pair<PropertyMap, LAS_property::Tag >` if the user wants to
   read a %LAS property as a scalar value `LAS_property::Tag::type` (for
   example, storing an `int` %LAS property into an `int` variable).

   - A `std::tuple<PropertyMap, Constructor,
   LAS_property::Tag...>` if the user wants to use one or several
   %LAS properties to construct a complex object (for example,
   storing 4 `unsigned short` %LAS properties into a %Color object
   that can for example be a `std::array<unsigned short,
   4>`). In that case, the second element of the tuple should be a
   functor that constructs the value type of `PropertyMap` from N
   objects of of type `LAS_property::Tag::type`.

   The %LAS standard defines a fixed set of properties accessible
   through the following tag classes:

   - `LAS_property::X` with type `double`
   - `LAS_property::Y` with type `double`
   - `LAS_property::Z` with type `double`
   - `LAS_property::Intensity` with type `unsigned short`
   - `LAS_property::Return_number` with type `unsigned char`
   - `LAS_property::Number_of_returns` with type `unsigned char`
   - `LAS_property::Scan_direction_flag` with type `unsigned char`
   - `LAS_property::Edge_of_flight_line` with type `unsigned char`
   - `LAS_property::Classification` with type `unsigned char`
   - `LAS_property::Synthetic_flag` with type `unsigned char`
   - `LAS_property::Keypoint_flag` with type `unsigned char`
   - `LAS_property::Withheld_flag` with type `unsigned char`
   - `LAS_property::Scan_angle` with type `double`
   - `LAS_property::User_data` with type `unsigned char`
   - `LAS_property::Point_source_ID` with type `unsigned short`
   - `LAS_property::Deleted_flag` with type `unsigned int`
   - `LAS_property::GPS_time` with type `double`
   - `LAS_property::R` with type `unsigned short`
   - `LAS_property::G` with type `unsigned short`
   - `LAS_property::B` with type `unsigned short`
   - `LAS_property::I` with type `unsigned short`

   \attention To read a binary file, the flag `std::ios::binary` must be set during the creation of the `ifstream`.

   \tparam OutputIteratorValueType type of objects that can be put in `PointOutputIterator`.
   It must be a model of `DefaultConstructible` and defaults to `value_type_traits<PointOutputIterator>::%type`.
   It can be omitted if the default is fine.
   \tparam PointOutputIterator iterator over output points.
   \tparam PropertyHandler handlers to recover properties.

   \returns `true` if reading was successful, `false` otherwise.

   \sa `make_las_point_reader()`

   \sa \ref IOStreamLAS
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename ... PropertyHandler>
bool read_LAS_with_properties(std::istream& is,
                              PointOutputIterator output,
                              PropertyHandler&& ... properties);



/**
   \ingroup PkgStreamSupportIoFuncsLAS

   \brief reads points (position only) using the \ref IOStreamLAS.

   Potential additional properties are ignored.

   \attention To read a binary file, the flag `std::ios::binary` must be set during the creation of the `ifstream`.

   \tparam OutputIteratorValueType type of objects that can be put in `PointOutputIterator`.
   It must be a model of `DefaultConstructible` and defaults to `value_type_traits<PointOutputIterator>::%type`.
   It can be omitted when the default is fine.
   \tparam PointOutputIterator iterator over output points.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param is input stream
   \param output output iterator over points
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if reading was successful, `false` otherwise.

   \sa `read_LAS_with_properties()`
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_LAS(std::istream& is,
              PointOutputIterator output,
              const CGAL_NP_CLASS& np = parameters::default_values());




/**
   \ingroup PkgStreamSupportIoFuncsLAS

   \brief reads points (position only) using the \ref IOStreamLAS.

   Potential additional properties are ignored.

   \tparam OutputIteratorValueType type of objects that can be put in `PointOutputIterator`.
   It must be a model of `DefaultConstructible` and defaults to `value_type_traits<PointOutputIterator>::%type`.
   It can be omitted when the default is fine.
   \tparam PointOutputIterator iterator over output points.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param filename name of the input file
   \param output output iterator over points
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if reading was successful, `false` otherwise.

   \sa `read_LAS_with_properties()`
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_LAS(const std::string& filename,
              PointOutputIterator output,
              const CGAL_NP_CLASS& np = parameters::default_values());


/**
 \ingroup PkgStreamSupportIoFuncsLAS

 \brief generates a %LAS property handler to write 3D points.

 \tparam PointMap the property map used to store points.

 \sa `write_LAS()`
 \sa \ref IOStreamLAS
*/
template <typename PointMap>
std::tuple<PointMap, LAS_property::X, LAS_property::Y, LAS_property::Z >
make_las_point_writer(PointMap point_map);

/**
   \ingroup PkgStreamSupportIoFuncsLAS

   \brief writes the range of `points` with properties to a .las stream.

   Properties are handled through a variadic list of property
   handlers. A `PropertyHandle` is a `std::pair<PropertyMap,
   LAS_property::Tag >` used to write a scalar value
   `LAS_property::Tag::type` as a %LAS property (for example,
   writing an `int` variable as an `int` %LAS property). An exception
   is used for points that are written using a `std::tuple` object.

   See documentation of `read_LAS_with_properties()` for the
   list of available `LAS_property::Tag` classes.

   \attention To write to a binary file, the flag `std::ios::binary` must be set during the creation of the `ofstream`.

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
                               PropertyHandler&& ... properties);

/**
   \ingroup PkgStreamSupportIoFuncsLAS

   \brief writes the range of `points` (positions only), using the \ref IOStreamLAS.

  \attention To write to a binary file, the flag `std::ios::binary` must be set during the creation of the `ofstream`.

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
template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_LAS(std::ostream& os,
               const PointRange& points,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PointRange>::value>* = nullptr
#endif
               );


/**
   \ingroup PkgStreamSupportIoFuncsLAS

   writes the range of `points` (positions only), using the \ref IOStreamLAS.

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
template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_LAS(const std::string& filename,
               const PointRange& points,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PointRange>::value>* = nullptr
#endif
               );

} // namespace IO

} // namespace CGAL


#ifdef CGAL_LINKED_WITH_LASLIB
#include <CGAL/IO/LAS/read_las_points.h>
#include <CGAL/IO/LAS/write_las_points.h>
#endif

#endif // CGAL_IO_LAS_H
