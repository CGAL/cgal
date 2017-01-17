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

#ifndef CGAL_READ_LAS_POINTS_H
#define CGAL_READ_LAS_POINTS_H

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Kernel_traits.h>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

#include <iostream>
#include <sstream>
#include <string>


namespace CGAL {


// PLY types:
// name        type        number of bytes
// ---------------------------------------
// char       character                 1
// uchar      unsigned character        1
// short      short integer             2
// ushort     unsigned short integer    2
// int        integer                   4
// uint       unsigned integer          4
// float      single-precision float    4
// double     double-precision float    8

/// \cond SKIP_IN_MANUAL

namespace internal {

  
} // namespace internal
  

/// \endcond



//===================================================================================
/// \ingroup PkgPointSetProcessing
/// Reads points (positions + normals, if available) from a .las
/// stream.
/// Potential additional point properties are ignored.
///
/// @tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
///         It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
/// @tparam OutputIterator iterator over output points.
/// @tparam PointPMap is a model of `WritablePropertyMap` with  value type `Point_3<Kernel>`.
///        It can be omitted if the value type of `OutputIterator` is convertible to `Point_3<Kernel>`.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with value type `Vector_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointPMap`.
///
/// @return true on success.

// This variant requires all parameters.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointPMap,
           typename NormalPMap,
           typename Kernel >
bool read_las_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_pmap, ///< property map: value_type of OutputIterator -> Vector_3.
                                 const Kernel& kernel) ///< geometric traits.
{

  return true;
}

/// @cond SKIP_IN_MANUAL
template < typename OutputIterator,
           typename PointPMap,
           typename NormalPMap,
           typename Kernel >
bool read_las_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_pmap, ///< property map: value_type of OutputIterator -> Vector_3.
                                 const Kernel& kernel) ///< geometric traits.
{

  return read_las_points_and_normals
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output,
                                                       point_pmap,
                                                       normal_pmap,
                                                       kernel);
}
//-----------------------------------------------------------------------------------
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointPMap,
           typename NormalPMap >
bool read_las_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return read_las_points_and_normals
    <OutputIteratorValueType>(stream,
                              output,
                              point_pmap,
                              normal_pmap,
                              Kernel());
}

template < typename OutputIterator,
           typename PointPMap,
           typename NormalPMap >
bool read_las_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{
  // just deduce value_type of OutputIterator
  return read_las_points_and_normals
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output,
                                                       point_pmap,
                                                       normal_pmap);
}
//-----------------------------------------------------------------------------------
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename NormalPMap >
bool read_las_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_las_points_and_normals
    <OutputIteratorValueType>(stream,
                              output,
                              make_identity_property_map(OutputIteratorValueType()),
                              normal_pmap);
}

template < typename OutputIterator,
           typename NormalPMap >
bool read_las_points_and_normals(std::istream& stream, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{
  // just deduce value_type of OutputIterator
  return read_las_points_and_normals
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output,
                                                       normal_pmap);
}
//-----------------------------------------------------------------------------------
/// @endcond


//===================================================================================
/// \ingroup PkgPointSetProcessing
/// Reads points (position only) from a .las stream.
/// Potential additional point properties (including normals).
///
/// @tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
///         It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
/// @tparam OutputIterator iterator over output points.
/// @tparam PointPMap is a model of `WritablePropertyMap` with  value_type `Point_3<Kernel>`.
///        It can be omitted if the value type of `OutputIterator` is convertible to `Point_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from  the value type of `PointPMap`.
///
/// @return `true` on success.

// This variant requires all parameters.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointPMap,
           typename Kernel >
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                     const Kernel& kernel) ///< geometric traits.
{
  // Calls read_las_points_and_normals() with a normal property map = boost::dummy_property_map
  return read_las_points_and_normals
    <OutputIteratorValueType>(stream,
                              output,
                              point_pmap,
                              boost::dummy_property_map(),
                              kernel);
}

/// @cond SKIP_IN_MANUAL
template < typename OutputIterator,
           typename PointPMap,
           typename Kernel >
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
                     const Kernel& kernel) ///< geometric traits.
{
  // just deduce value_type of OutputIterator
  return read_las_points
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output,
                                                       point_pmap,
                                                       kernel);
}
//-----------------------------------------------------------------------------------
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator,
           typename PointPMap >
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_pmap) ///< property map: value_type of OutputIterator -> Point_3.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return read_las_points
    <OutputIteratorValueType>(stream,
                              output,
                              point_pmap,
                              Kernel());
}

template < typename OutputIterator,
           typename PointPMap >
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_pmap) ///< property map: value_type of OutputIterator -> Point_3.
{
  // just deduce value_type of OutputIterator
  return read_las_points
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output,
                                                       point_pmap);
}
//-----------------------------------------------------------------------------------
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
//-----------------------------------------------------------------------------------
template < typename OutputIteratorValueType,
           typename OutputIterator >
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output) ///< output iterator over points.
{
  return read_las_points
    <OutputIteratorValueType>(stream,
                              output,
                              make_identity_property_map(OutputIteratorValueType())
                              );
}

template < typename OutputIterator>
bool read_las_points(std::istream& stream, ///< input stream.
                     OutputIterator output) ///< output iterator over points.
{
  // just deduce value_type of OutputIterator
  return read_las_points
    <typename value_type_traits<OutputIterator>::type>(stream,
                                                       output);
}
//-----------------------------------------------------------------------------------

/// @endcond


} //namespace CGAL

#endif // CGAL_READ_LAS_POINTS_H
