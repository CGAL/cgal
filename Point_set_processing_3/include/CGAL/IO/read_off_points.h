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

#ifndef CGAL_POINT_SET_PROCESSING_READ_OFF_POINTS_H
#define CGAL_POINT_SET_PROCESSING_READ_OFF_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/IO/io.h>
#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/Origin.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/iterator.h>
#include <CGAL/type_traits/is_iterator.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>

namespace CGAL {

namespace IO {

/**
   \ingroup PkgPointSetProcessing3IOOff

   \brief reads points (positions + normals, if available), using the \ref IOStreamOFF.

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

   \sa \ref IOStreamOFF
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(std::istream& is,
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
  typedef typename NP_helper::Geom_traits Kernel;
  typedef typename Kernel::FT FT;

  PointMap point_map = NP_helper::get_point_map(np);
  NormalMap normal_map = NP_helper::get_normal_map(np);

  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  // typedef typename value_type_traits<OutputIterator>::type Enriched_point;
  typedef OutputIteratorValueType Enriched_point;

  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;

  if(!is)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // scan points
  std::string signature;
  long pointsCount = 0, facesCount = 0, edgesCount = 0; // number of items in file
  int pointsRead = 0; // current number of points read
  int lineNumber = 0; // current line number
  std::string line;
  std::istringstream iss;
  while(getline(is,line))
  {
    iss.clear();
    iss.str(line);

    // Ignore empty lines and comments
    if (line.empty () || line[0] == '#')
      continue;

    ++lineNumber;

    // Reads file signature on first line
    if (lineNumber == 1)
    {
      if ( !(iss >> signature) || (signature != "OFF" && signature != "NOFF") )
      {
        // if wrong file format
        std::cerr << "Error line " << lineNumber << " of file (unexpected header)" << std::endl;
        return false;
      }
    }
    // Reads number of points on 2nd line
    else if (lineNumber == 2)
    {
      if ( !(iss >> pointsCount >> facesCount >> edgesCount) )
      {
        std::cerr << "Error line " << lineNumber << " of file (incorrect header format)" << std::endl;
        return false;
      }
    }

    // Reads 3D points on next lines
    else if (pointsRead < pointsCount)
    {
      // Reads position + normal...
      double x,y,z;
      double nx,ny,nz;
      if (iss >> IO::iformat(x) >> IO::iformat(y) >> IO::iformat(z))
      {
        Point point{FT(x), FT(y), FT(z)};
        Vector normal = CGAL::NULL_VECTOR;
        // ... + normal...
        if (iss >> IO::iformat(nx))
        {
          // In case we could read one number, we expect that there are two more
          if(iss  >> IO::iformat(ny) >> IO::iformat(nz)){
            normal = Vector(FT(nx),FT(ny),FT(nz));
          } else {
            std::cerr << "Error line " << lineNumber << " of file (incomplete normal coordinates)" << std::endl;
            return false;
          }
        }
        else if (signature == "NOFF")
        {
          std::cerr << "Error line " << lineNumber << " of file (expected normal coordinates)" << std::endl;
          return false;
        }

        Enriched_point pwn;
        put(point_map,  pwn, point);  // point_map[&pwn] = point
        put(normal_map, pwn, normal); // normal_map[&pwn] = normal

        *output++ = pwn;
        ++pointsRead;
      }
    }
  }

  if(is.eof())
    is.clear(is.rdstate() & ~std::ios_base::failbit); // set by getline

  return true;
}

/**
   \ingroup PkgPointSetProcessing3IOOff

   \brief reads points (positions + normals, if available), using the \ref IOStreamOFF.

   \tparam OutputIteratorValueType type of objects that can be put in `PointOutputIterator`.
   It must be a model of `DefaultConstructible` and defaults to `value_type_traits<PointOutputIterator>::%type`.
   It can be omitted when the default is fine.
   \tparam PointOutputIterator iterator over output points.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param fname input file name
   \param output output iterator over points
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

   \sa \ref IOStreamOFF
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const std::string& fname,
              PointOutputIterator output,
              const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
              , std::enable_if_t<CGAL::is_iterator<PointOutputIterator>::value>* = nullptr
#endif
              )
{
  std::ifstream is(fname);
  return read_OFF<OutputIteratorValueType>(is, output, np);
}

/// \cond SKIP_IN_MANUAL

// variants with default output iterator value type
template <typename OutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(std::istream& is, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values(),
              std::enable_if_t<CGAL::is_iterator<OutputIterator>::value>* = nullptr)
{
  return read_OFF<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

template <typename OutputIterator,typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const std::string& fname, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values(),
              std::enable_if_t<CGAL::is_iterator<OutputIterator>::value>* = nullptr)
{
  std::ifstream is(fname);
  return read_OFF<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

/// \endcond

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/// \cond SKIP_IN_MANUAL

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_off_points_and_normals(), please update your code")
bool read_off_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_map,  ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_map, ///< property map: value_type of OutputIterator -> Vector_3.
                                 const Kernel& /*kernel*/) ///< geometric traits.
{
  return IO::read_OFF<OutputIteratorValueType>(is, output,
                                               parameters::point_map(point_map)
                                                          .normal_map(normal_map)
                                                          .geom_traits(Kernel()));
}

template <typename OutputIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_off_points_and_normals(), please update your code")
bool read_off_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_map, ///< property map: value_type of OutputIterator -> Vector_3.
                                 const Kernel& kernel) ///< geometric traits.
{
  return IO::read_OFF<typename value_type_traits<OutputIterator>::type>(is, output,
                                                                        parameters::point_map(point_map)
                                                                                   .normal_map(normal_map)
                                                                                   .geom_traits(kernel));
}

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename NormalPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_off_points_and_normals(), please update your code")
bool read_off_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_OFF<OutputIteratorValueType>(is, output, parameters::point_map(point_map)
                                                                      .normal_map(normal_map));
}

template <typename OutputIterator,
          typename PointPMap,
          typename NormalPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_off_points_and_normals(), please update your code")
bool read_off_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                 NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_OFF<typename value_type_traits<OutputIterator>::type>(is, output,
                                                                        parameters::point_map(point_map)
                                                                                   .normal_map(normal_map));
}

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename NormalPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_off_points_and_normals(), please update your code")
bool read_off_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_OFF<OutputIteratorValueType>(is, output, parameters::normal_map(normal_map));
}

template <typename OutputIterator,
          typename NormalPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_off_points_and_normals(), please update your code")
bool read_off_points_and_normals(std::istream& is, ///< input stream.
                                 OutputIterator output, ///< output iterator over points.
                                 NormalPMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  return IO::read_OFF<typename value_type_traits<OutputIterator>::type>(is, output, parameters::normal_map(normal_map));
}

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_off_points(), please update your code")
bool read_off_points(std::istream& is, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                     const Kernel& kernel) ///< geometric traits.
{
  return read_off_points<OutputIteratorValueType>(is, output, parameters::point_map(point_map)
                                                                             .geom_traits(kernel));
}

template <typename OutputIterator,
          typename PointPMap,
          typename Kernel>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_off_points(), please update your code")
bool read_off_points(std::istream& is, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                     const Kernel& kernel) ///< geometric traits.
{
  return read_off_points<typename value_type_traits<OutputIterator>::type>(is, output,
                                                                           parameters::point_map(point_map)
                                                                                      .geom_traits (kernel));
}

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename PointPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_off_points(), please update your code")
bool read_off_points(std::istream& is, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_off_points<OutputIteratorValueType>(is, output, parameters::point_map (point_map));
}

template <typename OutputIterator, typename PointPMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::read_off_points(), please update your code")
bool read_off_points(std::istream& is, ///< input stream.
                     OutputIterator output, ///< output iterator over points.
                     PointPMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  return read_off_points<typename value_type_traits<OutputIterator>::type>(is, output, parameters::point_map(point_map));
}

/// \endcond

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool read_off_points(std::istream& is,
                                     OutputIterator output,
                                     const CGAL_NP_CLASS& np = parameters::default_values())
{
  return IO::read_OFF(is, output, np);
}

/// \cond SKIP_IN_MANUAL

template <typename OutputIteratorValueType, typename OutputIterator>
CGAL_DEPRECATED bool read_off_points(std::istream& is, OutputIterator output)
{
  return IO::read_OFF(is, output);
}

// variant with default output iterator value type
template <typename OutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool read_off_points(std::istream& is, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values())
{
  return IO::read_OFF(is, output, np);
}

/// \endcond

#endif //CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_POINT_SET_PROCESSING_READ_OFF_POINTS_H
