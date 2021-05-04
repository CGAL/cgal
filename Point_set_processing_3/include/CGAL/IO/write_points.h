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

#ifndef CGAL_POINT_SET_PROCESSING_WRITE_POINTS_H
#define CGAL_POINT_SET_PROCESSING_WRITE_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/IO/helpers.h>

#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/IO/write_xyz_points.h>
#ifdef CGAL_LINKED_WITH_LASLIB
#include <CGAL/IO/write_las_points.h>
#endif

#include <iostream>
#include <fstream>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#endif

namespace CGAL {

namespace IO {

/**
  \ingroup PkgPointSetProcessing3IO

  \brief writes the range of `points` with properties to a file.

  Supported file formats are the following:
  - \ref IOStreamOFF (`.off`)
  - \ref IOStreamPLY (`.ply`)
  - \ref IOStreamLAS (`.las`)
  - \ref IOStreamXYZ (`.xyz`)

  The format is detected from the filename extension (letter case is not important).

  \tparam PointRange is a model of `ConstRange`. The value type of
  its iterator is the key type of the named parameter `point_map`.
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the name of the output file.
  \param points the range of points that will be written.
  \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.

  \cgalNamedParamsBegin
    \cgalParamNBegin{point_map}
      \cgalParamDescription{a property map associating points to the elements of the point range}
      \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
      \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
    \cgalParamNEnd

    \cgalParamNBegin{normal_map}
      \cgalParamDescription{a property map associating normals to the elements of the point range}
      \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Vector_3`}
      \cgalParamDefault{If this parameter is omitted, normals are not written in the output file.}
    \cgalParamNEnd

    \cgalParamNBegin{geom_traits}
      \cgalParamDescription{an instance of a geometric traits class}
      \cgalParamType{a model of `Kernel`}
      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
    \cgalParamNEnd

     \cgalParamNBegin{use_binary_mode}
       \cgalParamDescription{indicates whether data should be written in binary (`true`) or in ASCII (`false`)}
       \cgalParamType{Boolean}
       \cgalParamDefault{`true`}
       \cgalParamExtra{This parameter is only relevant for `PLY` writing: the `OFF` and `XYZ` formats
                       are always ASCII, and the `LAS` format is always binary.}
     \cgalParamNEnd

    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`6`}
      \cgalParamExtra{This parameter is only meaningful while using ASCII encoding.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if writing was successful, `false` otherwise.
*/
template <typename PointRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_points(const std::string& fname,
                  const PointRange& points,
                  const CGAL_BGL_NP_CLASS& np,
#ifndef DOXYGEN_RUNNING
                  typename boost::enable_if<internal::is_Range<PointRange> >::type* = nullptr
#endif
                  )
{
  const std::string ext = internal::get_file_extension(fname);

  if(ext == "xyz" || ext == "pwn")
    return write_XYZ(fname, points, np);
  else if(ext == "off")
    return write_OFF(fname, points, np);
  else if(ext == "ply")
    return write_PLY(fname, points, np);
#ifdef CGAL_LINKED_WITH_LASLIB
  else if(ext == "las")
    return write_LAS(fname, points, np);

#endif
  return false;
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange>
bool write_points(const std::string& fname,const PointRange& points,
                  typename boost::enable_if<internal::is_Range<PointRange> >::type* = nullptr)
{
  return write_points(fname, points, parameters::all_default());
}

/// \endcond

} } // namespace CGAL::IO

#endif // CGAL_POINT_SET_PROCESSING_WRITE_POINTS_H
