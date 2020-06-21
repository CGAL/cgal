// Copyright (c) 2020 Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Maxime Gimeno

#ifndef CGAL_READ_POINTS_H
#define CGAL_READ_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/IO/polygon_soup_io.h>

#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/read_xyz_points.h>

#ifdef CGAL_LINKED_WITH_LASLIB
#include <CGAL/IO/read_las_points.h>
#endif

#include <fstream>
#include <string>

namespace CGAL {

/**
  \ingroup PkgPointSetProcessing3IO

  Reads the point set from an input file that can be either:

  - \link IOStreamXYZ XYZ \endlink
  - \link IOStreamOFF OFF \endlink
  - \link IOStreamPLY PLY \endlink
  - \link IOStreamLAS LAS \endlink

  \tparam OutputIteratorValueType type of objects that can be put in `PointOutputIterator`.
  It is default to `value_type_traits<PointOutputIterator>::%type` and can be omitted when the default is fine.
  \tparam PointOutputIterator iterator over output points.
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the name of the input file.
  \param output output iterator over points.
  \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

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

  \return `true` on success.
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename NamedParameters>
bool read_points(const std::string& fname,
                 PointOutputIterator output,
                 const NamedParameters& np)
{
  const std::string ext = IO::internal::get_file_extension(fname);

  if(ext == "xyz")
    return read_XYZ<OutputIteratorValueType>(fname, output, np);
  else if(ext == "off")
    return read_OFF<OutputIteratorValueType>(fname, output, np);
  else if(ext == "ply")
    return read_PLY<OutputIteratorValueType>(fname, output, np);
#ifdef CGAL_LINKED_WITH_LASLIB
  else if(ext == "las")
    return read_LAS<OutputIteratorValueType>(fname, output, np);
#endif

  return false;
}

// variant with default OutputIteratorType
template <typename OutputIterator, typename NamedParameters>
bool read_points(const std::string& fname, OutputIterator output, const NamedParameters& np)
{
  return read_points<typename value_type_traits<OutputIterator>::type>(fname, output, np);
}

// variant with default np
template <typename OutputIteratorValueType, typename OutputIterator>
bool read_points(const std::string& fname, OutputIterator output)
{
  return read_points<OutputIteratorValueType>(fname, output, parameters::all_default());
}

// variant with all default
template<typename OutputIterator>
bool read_points(const std::string& fname, OutputIterator output)
{
  return read_points<typename value_type_traits<OutputIterator>::type>(fname, output, parameters::all_default());
}

// variants with char*

template <typename OutputIteratorValueType, typename OutputIterator, typename NamedParameters>
bool read_points(const char* fname, OutputIterator output, const NamedParameters& np)
{
  return read_points<OutputIteratorValueType>(std::string(fname), output, np);
}

template <typename OutputIterator, typename NamedParameters>
bool read_points(const char* fname, OutputIterator output, const NamedParameters& np)
{
  return read_points<typename value_type_traits<OutputIterator>::type>(fname, output, np);
}

template <typename OutputIteratorValueType, typename OutputIterator>
bool read_points(const char* fname, OutputIterator output)
{
  return read_points<OutputIteratorValueType>(fname, output, parameters::all_default());
}

template<typename OutputIterator>
bool read_points(const char* fname, OutputIterator output)
{
  return read_points<typename value_type_traits<OutputIterator>::type>(fname, output, parameters::all_default());
}

} // namespace CGAL

#endif // CGAL_READ_POINTS_H
