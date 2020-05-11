// Copyright (c) 2020  Geometry Factory
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

#ifdef CGAL_LINKED_WITH_LASLIB
#include <CGAL/IO/read_las_points.h>
#endif
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/read_xyz_points.h>

namespace CGAL {

/**
  \ingroup PkgPointSetProcessing3IO
  Reads the point set from an input file that can be either:

  - XYZ
  - OFF
  - PLY
  - LAS
  \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
  It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
  \tparam OutputIterator iterator over output points.

  \param fname the name of the input file.
  \param output output iterator over points.
  \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

  \cgalNamedParamsBegin
    \cgalParamBegin{point_map} a model of `WritablePropertyMap` with value type `geom_traits::Point_3`.
    If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
    \cgalParamBegin{normal_map} a model of `ReadWritePropertyMap` with value type
    `geom_traits::Vector_3`. If this parameter is omitted, normals in the input stream are
    ignored.\cgalParamEnd
    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
  \cgalNamedParamsEnd

   \return true on success.
*/
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename NamedParameters>
bool read_points(const std::string& fname,
                    OutputIterator output,
                    const NamedParameters& np)
{
  if (fname.find(".xyz") != std::string::npos) {
    return read_XYZ<OutputIteratorValueType>(fname, output, np);
  }

  if (fname.find(".off") != std::string::npos) {
    return read_OFF<OutputIteratorValueType>(fname, output, np);
  }

  if (fname.find(".ply") != std::string::npos) {
    return read_PLY<OutputIteratorValueType>(fname, output, np);
  }

#ifdef CGAL_LINKED_WITH_LASLIB
  if (fname.find(".las") != std::string::npos) {
    return read_LAS<OutputIteratorValueType>(fname, output, np);
  }
#endif
  return false;
}

//variant with default OutputIteratorType
template <typename OutputIterator,
          typename NamedParameters>
bool read_points(const std::string& fname,
                    OutputIterator output,
                    const NamedParameters& np)
{
  return read_points<typename value_type_traits<OutputIterator>::type>(fname, output, np);
}

//variant with default np
template <typename OutputIteratorValueType,
          typename OutputIterator>
bool read_points(const std::string& fname,
                    OutputIterator output)
{
  return read_points<OutputIteratorValueType>(fname, output, parameters::all_default());
}

//variant with all default
template<typename OutputIterator>
bool read_points(const std::string& fname,
                 OutputIterator output)
{
  return read_points<typename value_type_traits<OutputIterator>::type>(fname, output, parameters::all_default());
}


//variants with char*

template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename NamedParameters>
bool read_points(const char* fname,
                    OutputIterator output,
                    const NamedParameters& np)
{
  return read_points<OutputIteratorValueType>(std::string(fname), output, np);
}

template <typename OutputIterator,
          typename NamedParameters>
bool read_points(const char* fname,
                    OutputIterator output,
                    const NamedParameters& np)
{
  return read_points<typename value_type_traits<OutputIterator>::type>(fname, output, np);
}

template <typename OutputIteratorValueType,
          typename OutputIterator>
bool read_points(const char* fname,
                    OutputIterator output)
{
  return read_points<OutputIteratorValueType>(fname, output, parameters::all_default());
}

template<typename OutputIterator>
bool read_points(const char* fname,
                 OutputIterator output)
{
  return read_points<typename value_type_traits<OutputIterator>::type>(fname, output, parameters::all_default());
}
}//end CGAL

#endif // CGAL_READ_POINTS_H
