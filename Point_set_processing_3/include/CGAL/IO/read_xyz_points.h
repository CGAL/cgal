// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_READ_XYZ_POINTS_H
#define CGAL_READ_XYZ_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>


#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Origin.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iostream>
#include <sstream>
#include <string>

namespace CGAL {

/**
   \ingroup PkgPointSetProcessing3IO
   Reads points (positions + normals, if available) from a .xyz ASCII stream.
   The function expects for each point a line with the x y z position,
   optionally followed by the nx ny nz normal.
   The first line may contain the number of points in the file.
   Empty lines and comments starting by # character are allowed.

   \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
   It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
   \tparam OutputIterator iterator over output points.

   \param stream input stream.
   \param output output iterator over points.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the poing range}
       \cgalParamType{a model of `ReadWritePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals in the input stream are ignored.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \return true on success.
*/
template <typename OutputIteratorValueType,
          typename OutputIterator,
#ifdef DOXYGEN_RUNNING
          typename NamedParameters
#else
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS
#endif
>
bool
read_xyz_points(
  std::istream& stream,
  OutputIterator output,
#ifdef DOXYGEN_RUNNING
    const NamedParameters& np)
#else
    const CGAL_BGL_NP_CLASS& np)
#endif
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3::Fake_point_range<OutputIteratorValueType> PointRange;

  // basic geometric types
  typedef typename CGAL::GetPointMap<PointRange, CGAL_BGL_NP_CLASS>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, CGAL_BGL_NP_CLASS>::type NormalMap;
  typedef typename Point_set_processing_3::GetK<PointRange, CGAL_BGL_NP_CLASS>::Kernel Kernel;

  bool has_normals = !(boost::is_same<NormalMap,
                       typename Point_set_processing_3::GetNormalMap<PointRange, CGAL_BGL_NP_CLASS>::NoMap>::value);

  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  NormalMap normal_map = choose_parameter<NormalMap>(get_parameter(np, internal_np::normal_map));

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
          put(point_map,  pwn, point);  // point_map[pwn] = point

          if (has_normals)
            put(normal_map, pwn, normal); // normal_map[pwn] = normal

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

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename OutputIteratorValueType,
          typename OutputIterator>
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output) ///< output iterator over points.
{
  return read_xyz_points<OutputIteratorValueType>
    (stream, output, CGAL::parameters::all_default());
}

// variant with default output iterator value type
template <typename OutputIterator,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output,
  const CGAL_BGL_NP_CLASS& np)
{
  return read_xyz_points<typename value_type_traits<OutputIterator>::type>
    (stream, output, np);
}

// variant with default NP and output iterator value type
template <typename OutputIterator>
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output)
{
  return read_xyz_points<typename value_type_traits<OutputIterator>::type>
    (stream, output, CGAL::parameters::all_default());
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_map,  ///< property map: value_type of OutputIterator -> Point_3.
  NormalPMap normal_map, ///< property map: value_type of OutputIterator -> Vector_3.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  return read_xyz_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map).
     geom_traits (Kernel()));
}

// deprecated API
template <typename OutputIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
  NormalPMap normal_map, ///< property map: value_type of OutputIterator -> Vector_3.
  const Kernel& kernel) ///< geometric traits.
{
  return read_xyz_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map).
     geom_traits (kernel));
}

// deprecated API
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename NormalPMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
  NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_xyz_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map));
}

// deprecated API
template <typename OutputIterator,
          typename PointPMap,
          typename NormalPMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
  NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_xyz_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map));
}

// deprecated API
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename NormalPMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_xyz_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::normal_map (normal_map));
}

// deprecated API
template <typename OutputIterator,
          typename NormalPMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool
read_xyz_points_and_normals(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return read_xyz_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::normal_map (normal_map));
}

// deprecated API
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points(), please update your code")
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
  const Kernel& kernel) ///< geometric traits.
{
  return read_xyz_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::point_map (point_map).
     geom_traits (kernel));
}

// deprecated API
template <typename OutputIterator,
          typename PointPMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points(), please update your code")
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
  const Kernel& kernel) ///< geometric traits.
{
  return read_xyz_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::point_map (point_map).
     geom_traits (kernel));
}

// deprecated API
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points(), please update your code")
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_xyz_points<OutputIteratorValueType>
    (stream, output,
     CGAL::parameters::point_map (point_map));
}

// deprecated API
template <typename OutputIterator,
          typename PointPMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points(), please update your code")
bool
read_xyz_points(
  std::istream& stream, ///< input stream.
  OutputIterator output, ///< output iterator over points.
  PointPMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_xyz_points<typename value_type_traits<OutputIterator>::type>
    (stream, output,
     CGAL::parameters::point_map (point_map));
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond


} //namespace CGAL

#endif // CGAL_READ_XYZ_POINTS_H
