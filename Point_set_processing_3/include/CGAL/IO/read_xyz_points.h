// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
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
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_READ_XYZ_POINTS_H
#define CGAL_READ_XYZ_POINTS_H

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>

#include <boost/version.hpp>
#if BOOST_VERSION >= 104000
  #include <boost/property_map/property_map.hpp>
#else
  #include <boost/property_map.hpp>
#endif

#include <iostream>
#include <sstream>
#include <string>

namespace CGAL {


//===================================================================================
/// \ingroup PkgPointSetProcessing
/// Reads points (positions + normals, if available) from a .xyz ASCII stream.
/// The function expects for each point a line with the x y z position,
/// optionally followed by the nx ny nz normal.
/// The first line may contain the number of points in the file.
/// Empty lines and comments starting by # character are allowed.
///
/// @tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
///         It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
/// @tparam OutputIterator iterator over output points.
/// @tparam PointPMap is a model of `WritablePropertyMap` with  value type `Point_3<Kernel>`.
///        It can be omitted if the value type of OutputIterator value_type is convertible to `Point_3<Kernel>`.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with value type  `Vector_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointPMap`.
///
/// @return true on success.

// This variant requires all parameters.
//-----------------------------------------------------------------------------------
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel
>
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap,  ///< property map: value_type of OutputIterator -> Point_3.
  NormalPMap normal_pmap, ///< property map: value_type of OutputIterator -> Vector_3.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  //typedef typename value_type_traits<OutputIterator>::type Enriched_point;
  typedef OutputIteratorValueType Enriched_point;
  
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;

  if(!stream)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // scan points
  long pointsCount; // number of points in file
  int lineNumber = 0; // line counter
  std::string line; // line buffer
  std::istringstream iss;

  while(getline(stream,line))
  {
    // position + normal
    FT x,y,z;
    FT nx,ny,nz;

    lineNumber++;

    // Trims line buffer
    line.erase(line.find_last_not_of (" ")+1);
    line.erase(0, line.find_first_not_of (" "));

    // Skips comment or empty line...
    if (line.length() == 0 || line[0] == '#')
    {
      continue;
    }
    // ...or reads position...
    else {
      iss.clear();
      iss.str(line);
      if (iss >> iformat(x) >> iformat(y) >> iformat(z))
        {
          Point point(x,y,z);
          Vector normal = CGAL::NULL_VECTOR;
          // ... + normal...
          if (iss >> iformat(nx))
            {
              // In case we could read one number, we expect that there are two more
              if(iss  >> iformat(ny) >> iformat(nz)){
                normal = Vector(nx,ny,nz);
              } else {
                std::cerr << "Error line " << lineNumber << " of file" << std::endl;
                return false;
              }
            }
          Enriched_point pwn;
        #ifdef CGAL_USE_PROPERTY_MAPS_API_V1
          put(point_pmap,  &pwn, point);  // point_pmap[&pwn] = point
          put(normal_pmap, &pwn, normal); // normal_pmap[&pwn] = normal
        #else
          put(point_pmap,  pwn, point);  // point_pmap[pwn] = point
          put(normal_pmap, pwn, normal); // normal_pmap[pwn] = normal
        #endif
          *output++ = pwn;
          continue;
        } 
      
    }
    // ...or skips number of points on first line (optional)
    if (lineNumber == 1 && std::istringstream(line) >> pointsCount)
    {
      continue;
    }
    else // if wrong file format
    {
      std::cerr << "Error line " << lineNumber << " of file" << std::endl;
      return false;
    }
  }

  return true;
}

/// @cond SKIP_IN_MANUAL
template <typename OutputIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel
>
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
  NormalPMap normal_pmap, ///< property map: value_type of OutputIterator -> Vector_3.
  const Kernel& kernel) ///< geometric traits.
{
  // just deduce value_type of OutputIterator
  return read_xyz_points_and_normals
    <typename value_type_traits<OutputIterator>::type>(
    stream,
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
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename NormalPMap
>
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
  NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return read_xyz_points_and_normals
    <OutputIteratorValueType>(
    stream,
    output,
    point_pmap,
    normal_pmap,
    Kernel());
}

template <typename OutputIterator,
          typename PointPMap,
          typename NormalPMap
>
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
  NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{
  // just deduce value_type of OutputIterator
  return read_xyz_points_and_normals
    <typename value_type_traits<OutputIterator>::type>(
    stream,
    output,
    point_pmap,
    normal_pmap);
}
//-----------------------------------------------------------------------------------
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
//-----------------------------------------------------------------------------------
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename NormalPMap
>
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_xyz_points_and_normals
    <OutputIteratorValueType>(
    stream,
    output,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(output),
#else
    make_identity_property_map(OutputIteratorValueType()),
#endif
    normal_pmap);
}

template <typename OutputIterator,
          typename NormalPMap
>
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  NormalPMap normal_pmap) ///< property map: value_type of OutputIterator -> Vector_3.
{
  // just deduce value_type of OutputIterator
  return read_xyz_points_and_normals
    <typename value_type_traits<OutputIterator>::type>(
    stream,
    output,
    normal_pmap);
}
//-----------------------------------------------------------------------------------
/// @endcond


//===================================================================================
/// \ingroup PkgPointSetProcessing
/// Reads points (positions only) from a .xyz ASCII stream.
/// The function expects for each point a line with the x y z position.
/// If the position is followed by the nx ny nz normal, then the normal will be ignored.
/// The first line may contain the number of points in the file.
/// Empty lines and comments starting by # character are allowed.
///
/// @tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
///         It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
/// @tparam OutputIterator iterator over output points.
/// @tparam PointPMap is a model of `WritablePropertyMap` with value type  `Point_3<Kernel>`.
///        It can be omitted if the value type of `OutputIterator` is convertible to `Point_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointPMap`.
///
/// @return true on success.

// This variant requires all parameters.
//-----------------------------------------------------------------------------------
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename Kernel
>
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
  const Kernel& kernel) ///< geometric traits.
{
  // Calls read_xyz_points_and_normals() with a normal property map = boost::dummy_property_map
  return read_xyz_points_and_normals
    <OutputIteratorValueType>(
    stream,
    output,
    point_pmap,
    boost::dummy_property_map(),
    kernel);
}

/// @cond SKIP_IN_MANUAL
template <typename OutputIterator,
          typename PointPMap,
          typename Kernel
>
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map: value_type of OutputIterator -> Point_3.
  const Kernel& kernel) ///< geometric traits.
{
  // just deduce value_type of OutputIterator
  return read_xyz_points
    <typename value_type_traits<OutputIterator>::type>(
    stream,
    output,
    point_pmap,
    kernel);
}
/// @endcond
//-----------------------------------------------------------------------------------

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
//-----------------------------------------------------------------------------------
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap
>
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap) ///< property map: value_type of OutputIterator -> Point_3.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return read_xyz_points
    <OutputIteratorValueType>(
    stream,
    output,
    point_pmap,
    Kernel());
}

template <typename OutputIterator,
          typename PointPMap
>
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap) ///< property map: value_type of OutputIterator -> Point_3.
{
  // just deduce value_type of OutputIterator
  return read_xyz_points
    <typename value_type_traits<OutputIterator>::type>(
    stream,
    output,
    point_pmap);
}
//-----------------------------------------------------------------------------------
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
//-----------------------------------------------------------------------------------
template <typename OutputIteratorValueType,
          typename OutputIterator
>
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output) ///< output iterator over points.
{
  return read_xyz_points
    <OutputIteratorValueType>(
    stream,
    output,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(output)
#else
    make_identity_property_map(OutputIteratorValueType())
#endif
    );
}

template <typename OutputIterator
>
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output) ///< output iterator over points.
{
  // just deduce value_type of OutputIterator
  return read_xyz_points
    <typename value_type_traits<OutputIterator>::type>(
    stream,
    output);
}
//-----------------------------------------------------------------------------------
/// @endcond


} //namespace CGAL

#endif // CGAL_READ_XYZ_POINTS_H
