// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

#ifndef CGAL_READ_OFF_POINTS_H
#define CGAL_READ_OFF_POINTS_H

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
/// Reads points (positions + normals, if available) from a .off ASCII stream.
/// The function expects for each point a line with the x y z position,
/// optionally followed by the nx ny nz normal.
/// Faces are ignored.
///
/// @commentheading Template Parameters:
/// @param OutputIterator iterator over output points.
/// @param PointPMap is a model of boost::WritablePropertyMap with a value_type = Point_3<Kernel>.
///        It can be omitted if OutputIterator value_type is convertible to Point_3<Kernel>.
/// @param NormalPMap is a model of boost::WritablePropertyMap with a value_type = Vector_3<Kernel>.
/// @param Kernel Geometric traits class.
///        It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return true on success.

// This variant requires all parameters.
template <typename OutputIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel
>
bool
read_off_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map OutputIterator -> Point_3.
  NormalPMap normal_pmap, ///< property map OutputIterator -> Vector_3.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  typedef typename value_type_traits<OutputIterator>::type Enriched_point;

  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;

  if(!stream)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // scan points
  long pointsCount = 0, facesCount = 0, edgesCount = 0; // number of items in file
  int pointsRead = 0; // current number of points read
  int lineNumber = 0; // current line number
  std::string line;
  while(getline(stream,line))
  {
    lineNumber++;

    // Reads file signature on first line
    if (lineNumber == 1)
    {
      std::string signature;
      std::istringstream iss(line);
      if ( !(iss >> signature)
        || (signature != "OFF" && signature != "NOFF") )
      {
        // if wrong file format
        std::cerr << "Incorrect file format line " << lineNumber << " of file" << std::endl;
        return false;
      }
    }

    // Reads number of points on 2nd line
    else if (lineNumber == 2)
    {
      std::istringstream iss(line);
      if ( !(iss >> pointsCount >> facesCount >> edgesCount) )
      {
        std::cerr << "Error line " << lineNumber << " of file" << std::endl;
        return false;
      }
    }

    // Reads 3D points on next lines
    else if (pointsRead < pointsCount)
    {
      // Reads position + normal...
      double x,y,z;
      double nx,ny,nz;
      if (std::istringstream(line) >> x >> y >> z >> nx >> ny >> nz)
      {
        Point point(x,y,z);
        Vector normal(nx,ny,nz);
        Enriched_point pwn;
        put(point_pmap,  &pwn, point);  // point_pmap[&pwn] = point
        put(normal_pmap, &pwn, normal); // normal_pmap[&pwn] = normal
        *output++ = pwn;
        pointsRead++;
      }
      // ...or read only position...
      else if (std::istringstream(line) >> x >> y >> z)
      {
        Point point(x,y,z);
        Vector normal = CGAL::NULL_VECTOR;
        Enriched_point pwn;
        put(point_pmap,  &pwn, point);  // point_pmap[&pwn] = point
        put(normal_pmap, &pwn, normal); // normal_pmap[&pwn] = normal
        *output++ = pwn;
        pointsRead++;
      }
      // ...or skip comment line
    }
    // Skip remaining lines
  }

  return true;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
template <typename OutputIterator,
          typename PointPMap,
          typename NormalPMap
>
bool
read_off_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map OutputIterator -> Point_3.
  NormalPMap normal_pmap) ///< property map OutputIterator -> Vector_3.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return read_off_points_and_normals(
    stream,
    output,
    point_pmap,
    normal_pmap,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map.
template <typename OutputIterator,
          typename NormalPMap
>
bool
read_off_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  NormalPMap normal_pmap) ///< property map OutputIterator -> Vector_3.
{
  return read_off_points_and_normals(
    stream,
    output,
    make_dereference_property_map(output),
    normal_pmap);
}
/// @endcond


//===================================================================================
/// Reads points (position only) from a .off ASCII stream.
/// The function expects for each point a line with the x y z position.
/// If the position is followed by the nx ny nz normal, then the normal will be ignored.
/// Faces are ignored.
///
/// @commentheading Template Parameters:
/// @param OutputIterator iterator over output points.
/// @param PointPMap is a model of boost::WritablePropertyMap with a value_type = Point_3<Kernel>.
///        It can be omitted if OutputIterator value_type is convertible to Point_3<Kernel>.
/// @param Kernel Geometric traits class.
///        It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return true on success.

// This variant requires all parameters.
template <typename OutputIterator,
          typename PointPMap,
          typename Kernel
>
bool
read_off_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map OutputIterator -> Point_3.
  const Kernel& kernel) ///< geometric traits.
{
  // Calls read_off_points_and_normals() with a normal property map = boost::dummy_property_map
  return read_off_points_and_normals(
    stream,
    output,
    point_pmap,
    boost::dummy_property_map(),
    kernel);
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
template <typename OutputIterator,
          typename PointPMap
>
bool
read_off_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap) ///< property map OutputIterator -> Point_3.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return read_off_points(
    stream,
    output,
    point_pmap,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map.
template <typename OutputIterator
>
bool
read_off_points(
  std::istream& stream, ///< input stream.
  OutputIterator output) ///< output iterator over points.
{
  return read_off_points(
    stream,
    output,
    make_dereference_property_map(output));
}
/// @endcond


} //namespace CGAL

#endif // CGAL_READ_OFF_POINTS_H
