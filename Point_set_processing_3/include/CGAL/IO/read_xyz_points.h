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
#include <CGAL/Origin.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/type_traits/is_iterator.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

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
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(std::istream& is,
              OutputIterator output,
              const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3::Fake_point_range<OutputIteratorValueType> PointRange;

  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, CGAL_NP_CLASS> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;
  typedef typename NP_helper::Geom_traits Kernel;

  PointMap point_map = NP_helper::get_point_map(np);
  NormalMap normal_map = NP_helper::get_normal_map(np);

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
    else
    {
      iss.clear();
      iss.str(line);
      if (iss >> IO::iformat(x) >> IO::iformat(y) >> IO::iformat(z))
      {
        Point point(x,y,z);
        Vector normal = CGAL::NULL_VECTOR;
        // ... + normal...
        if (iss >> IO::iformat(nx))
        {
          // In case we could read one number, we expect that there are two more
          if(iss >> IO::iformat(ny) >> IO::iformat(nz)){
            normal = Vector(nx,ny,nz);
          } else {
            std::cerr << "Error line " << lineNumber << " of file (incomplete normal coordinates)" << std::endl;
            return false;
          }
        }

        Enriched_point pwn;
        put(point_map,  pwn, point);  // point_map[pwn] = point

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
      std::cerr << "Error line " << lineNumber << " of file (expected point coordinates)" << std::endl;
      return false;
    }
  }

  if(is.eof())
    is.clear(is.rdstate() & ~std::ios_base::failbit); // set by getline

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
           typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(const std::string& fname,
              OutputIterator output,
              const CGAL_NP_CLASS& np = parameters::default_values())
{
  std::ifstream is(fname);
  return read_XYZ<OutputIteratorValueType>(is, output, np);
}

/// \cond SKIP_IN_MANUAL

// variants with default output iterator value type
template <typename OutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(std::istream& is, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values(),
              std::enable_if_t<CGAL::is_iterator<OutputIterator>::value>* = nullptr)
{
  return read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

template <typename OutputIterator,typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(const std::string& fname, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values())
{
  std::ifstream is(fname);
  return read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output, np);
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
          typename CGAL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool read_xyz_points(std::istream& is,
                                     OutputIterator output,
                                     const CGAL_NP_CLASS& np = parameters::default_values())
{
  return IO::read_XYZ(is, output, np);
}

/// \cond SKIP_IN_MANUAL
template <typename OutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool read_xyz_points(std::istream& is,
                                     OutputIterator output,
                                     const CGAL_NP_CLASS& np = parameters::default_values())
{
  return IO::read_XYZ<typename value_type_traits<OutputIterator>::type>(is, output, np);
}
/// \endcond

#endif //CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_POINT_SET_PROCESSING_READ_XYZ_POINTS_H
