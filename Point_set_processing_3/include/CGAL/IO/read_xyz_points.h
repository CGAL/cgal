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

#ifndef CGAL_POINT_SET_PROCESSING_READ_XYZ_POINTS_H
#define CGAL_POINT_SET_PROCESSING_READ_XYZ_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Origin.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/is_iterator.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#define CGAL_DEPRECATED
#endif

namespace CGAL {

namespace IO {

/**
   \ingroup PkgPointSetProcessing3IOXyz

   \brief reads points (positions + normals, if available), using the \ref IOStreamXYZ.

   \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
   It must be a model of `DefaultConstructible` and defaults to `value_type_traits<OutputIterator>::%type`.
   It can be omitted when the default is fine.
   \tparam OutputIterator iterator over output points.
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

   \sa \ref IOStreamXYZ
*/
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(std::istream& is,
              OutputIterator output,
              const CGAL_BGL_NP_CLASS& np)
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

  if(!is)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // scan points
  long pointsCount; // number of points in file
  int lineNumber = 0; // line counter
  std::string line; // line buffer
  std::istringstream iss;

  while(getline(is,line))
  {
    // position + normal
    FT x,y,z;
    FT nx,ny,nz;

    ++lineNumber;

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

/**
   \ingroup PkgPointSetProcessing3IOXyz

   \brief reads points (positions + normals, if available), using the \ref IOStreamXYZ.

   \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
   It must be a model of `DefaultConstructible` and defaults to `value_type_traits<OutputIterator>::%type`.
   It can be omitted when the default is fine.
   \tparam OutputIterator iterator over output points.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param fname input file name.
   \param output output iterator over points.
   \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.

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

   \sa \ref IOStreamXYZ
*/
template <typename OutputIteratorValueType,
          typename OutputIterator,
           typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(const std::string& fname,
              OutputIterator output,
              const CGAL_BGL_NP_CLASS& np)
{
  std::ifstream is(fname);
  return read_XYZ<OutputIteratorValueType>(is, output, np);
}

/// \cond SKIP_IN_MANUAL

// variants with default NP
template <typename OutputIteratorValueType, typename OutputIterator>
bool read_XYZ(std::istream& is, OutputIterator output)
{
  return read_XYZ<OutputIteratorValueType>(is, output, parameters::all_default());
}

template <typename OutputIteratorValueType, typename OutputIterator>
bool read_XYZ(const std::string& fname, OutputIterator output)
{
  return read_XYZ<OutputIteratorValueType>(fname, output, parameters::all_default());
}

// variants with default output iterator value type
template <typename OutputIterator, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(std::istream& is, OutputIterator output, const CGAL_BGL_NP_CLASS& np)
{
  return read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

template <typename OutputIterator,typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(const std::string& fname, OutputIterator output, const CGAL_BGL_NP_CLASS& np)
{
  std::ifstream is(fname);
  return read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

// variants with default NP and output iterator value type
template <typename OutputIterator>
bool read_XYZ(std::istream& is,
              OutputIterator output,
              typename std::enable_if<CGAL::is_iterator<OutputIterator>::value>::type* = nullptr)
{
  return read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output, parameters::all_default());
}

template <typename OutputIterator>
bool read_XYZ(const std::string& fname, OutputIterator output)
{
  return read_XYZ<typename value_type_traits<OutputIterator>::type>(fname, output, parameters::all_default());
}

} // namespace IO

/// \endcond

#ifndef CGAL_NO_DEPRECATED_CODE

/// \cond SKIP_IN_MANUAL

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool read_xyz_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_map,  ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_map, ///< property map: value_type of OutputIterator -> Vector_3.
                                 const Kernel& /*kernel*/) ///< geometric traits.
{
  return IO::read_XYZ<OutputIteratorValueType>(is, output,
                                               parameters::point_map(point_map)
                                                          .normal_map(normal_map)
                                                          .geom_traits(Kernel()));
}

template <typename OutputIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool read_xyz_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_map, ///< property map: value_type of OutputIterator -> Vector_3.
                                 const Kernel& kernel) ///< geometric traits.
{
  return IO::read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output,
                                                                        parameters::point_map(point_map)
                                                                                   .normal_map(normal_map)
                                                                                   .geom_traits(kernel));
}

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename NormalPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool read_xyz_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_XYZ<OutputIteratorValueType>(is, output,
                                               parameters::point_map(point_map)
                                                          .normal_map(normal_map));
}

template <typename OutputIterator,
          typename PointPMap,
          typename NormalPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool read_xyz_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output,
                                                                        parameters::point_map(point_map)
                                                                                   .normal_map(normal_map));
}

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename NormalPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool read_xyz_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_XYZ<OutputIteratorValueType>(is, output, parameters::normal_map(normal_map));
}

template <typename OutputIterator,
          typename NormalPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points_and_normals(), please update your code")
bool read_xyz_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output,
                                                                        parameters::normal_map(normal_map));
}

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename Kernel>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points(), please update your code")
bool read_xyz_points(std::istream& is, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                     const Kernel& kernel) ///< geometric traits.
{
  return IO::read_XYZ<OutputIteratorValueType>(is, output,
                                               parameters::point_map(point_map)
                                                          .geom_traits(kernel));
}

template <typename OutputIterator,
          typename PointPMap,
          typename Kernel>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points(), please update your code")
bool read_xyz_points(std::istream& is, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                     const Kernel& kernel) ///< geometric traits.
{
  return IO::read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output,
                                                                        parameters::point_map(point_map)
                                                                                   .geom_traits(kernel));
}

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points(), please update your code")
bool read_xyz_points(std::istream& is, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return IO::read_XYZ<OutputIteratorValueType>(is, output, parameters::point_map(point_map));
}

template <typename OutputIterator,
          typename PointPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_xyz_points(), please update your code")
bool read_xyz_points(std::istream& is, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return IO::read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output,
                                                                        parameters::point_map(point_map));
}

/// \endcond

/**
   \ingroup PkgPointSetProcessing3IODeprecated

   \deprecated This function is deprecated since \cgal 5.3,
               \link PkgPointSetProcessing3IOXyz `CGAL::IO::read_XYZ()` \endlink should be used instead.

   \returns `true` if reading was successful, `false` otherwise.
*/
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool read_xyz_points(std::istream& is,
                                     OutputIterator output,
                                     const CGAL_BGL_NP_CLASS& np)
{
  return IO::read_XYZ(is, output, np);
}

/// \cond SKIP_IN_MANUAL

template <typename OutputIteratorValueType,
          typename OutputIterator>
CGAL_DEPRECATED bool read_xyz_points(std::istream& is,
                                     OutputIterator output)
{
  return IO::read_XYZ<OutputIteratorValueType>(is, output, parameters::all_default());
}

template <typename OutputIterator,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool read_xyz_points(std::istream& is,
                                     OutputIterator output,
                                     const CGAL_BGL_NP_CLASS& np)
{
  return IO::read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

template <typename OutputIterator>
CGAL_DEPRECATED bool read_xyz_points(std::istream& is,
                                     OutputIterator output)
{
  return IO::read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output, parameters::all_default());
}

/// \endcond

#endif //CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_POINT_SET_PROCESSING_READ_XYZ_POINTS_H
