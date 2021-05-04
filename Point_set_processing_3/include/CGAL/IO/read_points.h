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

#ifndef CGAL_POINT_SET_PROCESSING_READ_POINTS_H
#define CGAL_POINT_SET_PROCESSING_READ_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/IO/helpers.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/read_xyz_points.h>

#ifdef CGAL_LINKED_WITH_LASLIB
#include <CGAL/IO/read_las_points.h>
#endif

#include <fstream>
#include <string>

namespace CGAL {

namespace IO {

/**
  \ingroup PkgPointSetProcessing3IO

  \brief reads the point set from an input file.

  Supported file formats are the following:
  - \ref IOStreamOFF (`.off`)
  - \ref IOStreamPLY (`.ply`)
  - \ref IOStreamLAS (`.las`)
  - \ref IOStreamXYZ (`.xyz`)

  The format is detected from the filename extension (letter case is not important).

  \tparam OutputIteratorValueType type of objects that can be put in `PointOutputIterator`.
  It must be a model of `DefaultConstructible` and defaults to `value_type_traits<PointOutputIterator>::%type`.
  It can be omitted when the default is fine.
  \tparam PointOutputIterator iterator over output points.
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the name of the input file.
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

     \cgalParamNBegin{use_binary_mode}
       \cgalParamDescription{indicates whether data should be read in binary (`true`) or in ASCII (`false`)}
       \cgalParamType{Boolean}
       \cgalParamDefault{`true`}
       \cgalParamExtra{This parameter is only relevant for `PLY` reading: the `OFF` and `XYZ` formats
                       are always ASCII, and the `LAS` format is always binary.}
     \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if reading was successful, `false` otherwise.
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename NamedParameters>
bool read_points(const std::string& fname,
                 PointOutputIterator output,
                 const NamedParameters& np)
{
  const std::string ext = internal::get_file_extension(fname);

  if(ext == "xyz" || ext == "pwn")
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

/// \cond SKIP_IN_MANUAL

// variant with default OutputIteratorType
template <typename OutputIterator, typename NamedParameters>
bool read_points(const std::string& fname, OutputIterator output, const NamedParameters& np)
{
  return read_points<typename value_type_traits<OutputIterator>::type>(fname, output, np);
}

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

/// \endcond

} } // namespace CGAL::IO

#endif // CGAL_POINT_SET_PROCESSING_READ_POINTS_H
