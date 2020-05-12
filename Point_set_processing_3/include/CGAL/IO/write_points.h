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

#ifndef CGAL_WRITE_POINTS_H
#define CGAL_WRITE_POINTS_H


#ifdef CGAL_LINKED_WITH_LASLIB
#include <CGAL/IO/write_las_points.h>
#endif
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/IO/write_xyz_points.h>

namespace CGAL {
/**
  \ingroup PkgPointSetProcessing3IO
   Saves the range of `points` with properties to a
   file that can be either:
  - XYZ
  - OFF
  - PLY
  - LAS
   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.

  \param fname the name of the input file.
  \param points the range of points that will be written.
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
template <typename PointRange,
          #ifdef DOXYGEN_RUNNING
          typename NamedParameters
          #else
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS
          #endif
          >
bool write_points(const std::string& fname,
                  const PointRange& points,
                    #ifdef DOXYGEN_RUNNING
                    const NamedParameters& np)
                #else
                    const CGAL_BGL_NP_CLASS& np)
                #endif
                {
  if (fname.find(".xyz") != std::string::npos) {
    return write_XYZ(fname, points, np);
  }

  if (fname.find(".off") != std::string::npos) {
    return write_OFF(fname, points, np);
  }

  if (fname.find(".ply") != std::string::npos) {
    return write_PLY(fname, points, np);
  }

#ifdef CGAL_LINKED_WITH_LASLIB
  if (fname.find(".las") != std::string::npos) {
    return write_LAS(fname, points, np);
  }
#endif
  return false;
}

template <typename PointRange,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_points(const char* fname,
                  const PointRange& points,
                    const CGAL_BGL_NP_CLASS& np)
{
  return write_points(std::string(fname), points, np);
}

template <typename PointRange>
bool write_points(const std::string& fname,const PointRange& points)
{
  return write_points(fname, points, parameters::all_default());
}

template <typename PointRange>
bool write_points(const char* fname,const PointRange& points)
{
  return write_points(fname, points, parameters::all_default());
}

}//end CGAL
#endif // CGAL_WRITE_POINTS_H
