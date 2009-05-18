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

#ifndef CGAL_READ_XYZ_POINTS_H
#define CGAL_READ_XYZ_POINTS_H

#include <CGAL/point_set_property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>

#include <boost/property_map.hpp>

#include <iostream>
#include <sstream>
#include <string>

CGAL_BEGIN_NAMESPACE


//===================================================================================
/// Reads points (positions + normals, if available) from a .xyz ASCII stream.
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
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map OutputIterator -> Point_3.
  NormalPMap normal_pmap, ///< property map OutputIterator -> Vector_3.
  const Kernel& kernel) ///< geometric traits.
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
  long pointsCount; // number of points in file
  int lineNumber = 0;
  std::string line;
  while(getline(stream,line))
  {
    lineNumber++;

    // Reads position + normal...
    double x,y,z;
    double nx,ny,nz;
    if (std::istringstream(line) >> x >> y >> z >> nx >> ny >> nz)
    {
      Point point(x,y,z);
      Vector normal(nx,ny,nz);
      //*output++ = Enriched_point(point,normal);
      Enriched_point pwn;
      put(point_pmap,  &pwn, point);  // point_pmap[&pwn] = point
      put(normal_pmap, &pwn, normal); // normal_pmap[&pwn] = normal
      *output++ = pwn;
    }
    // ...or read only position...
    else if (std::istringstream(line) >> x >> y >> z)
    {
      Point point(x,y,z);
      //*output++ = point;
      Enriched_point pwn;
      put(point_pmap,  &pwn, point);  // point_pmap[&pwn] = point
      *output++ = pwn;
    }
    // ...or skip number of points on first line (optional)
    else if (lineNumber == 1 && std::istringstream(line) >> pointsCount)
    {
      continue;
    }
    else
    {
      std::cerr << "Error line " << lineNumber << " of file" << std::endl;
      return false;
    }
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
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map OutputIterator -> Point_3.
  NormalPMap normal_pmap) ///< property map OutputIterator -> Vector_3.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return read_xyz_points_and_normals(
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
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  NormalPMap normal_pmap) ///< property map OutputIterator -> Vector_3.
{
  return read_xyz_points_and_normals(
    stream,
    output,
    make_dereference_property_map(output),
    normal_pmap);
}
/// @endcond


//===================================================================================
/// Reads points (positions only) from a .xyz ASCII stream.
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
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap, ///< property map OutputIterator -> Point_3.
  const Kernel& kernel) ///< geometric traits.
{
  // Calls read_xyz_points_and_normals() with a normal property map = boost::dummy_property_map
  return read_xyz_points_and_normals(
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
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_pmap) ///< property map OutputIterator -> Point_3.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return read_xyz_points(
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
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output) ///< output iterator over points.
{
  return read_xyz_points(
    stream,
    output,
    make_dereference_property_map(output));
}
/// @endcond


CGAL_END_NAMESPACE

#endif // CGAL_READ_XYZ_POINTS_H
