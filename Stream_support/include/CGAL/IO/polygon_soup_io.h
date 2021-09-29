// Copyright (c) 2020  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_IO_READ_POLYGON_SOUP_H
#define CGAL_IO_READ_POLYGON_SOUP_H

#include <CGAL/IO/3MF.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/IO/STL.h>
#include <CGAL/IO/VTK.h>
#include <CGAL/IO/GOCAD.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {

namespace IO {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

/*!
 * \ingroup IOstreamFunctions
 *
 * \brief reads a polygon soup from a file.
 *
 * Supported file formats are the following:
 * - \ref IOStreamOFF (`.off`)
 * - \ref IOStreamOBJ (`.obj`)
 * - \ref IOStreamSTL (`.stl`)
 * - \ref IOStreamPLY (`.ply`)
 * - \ref IOStreamGocad (`.ts`)
 * - \ref IOStreamVTK (`.vtp`)
 *
 * The format is detected from the filename extension (letter case is not important).
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
 * \tparam PolygonRange a model of the concepts `SequenceContainer` and `BackInsertionSequence`
 *                      whose `value_type` is itself a model of the concepts `SequenceContainer`
 *                      and `BackInsertionSequence` whose `value_type` is an unsigned integer type
 *                      convertible to `std::size_t`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the name of the file.
 * \param points points of the soup of polygons
 * \param polygons each element in the range describes a polygon using the indices of the vertices.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{indicates whether output warnings and error messages should be printed or not.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if reading was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_polygon_soup(const std::string& fname,
                       PointRange& points,
                       PolygonRange& polygons,
                       const CGAL_BGL_NP_CLASS& np)
{
  const bool verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false);

  const std::string ext = internal::get_file_extension(fname);
  if(ext == std::string())
  {
    if(verbose)
      std::cerr << "Error: cannot read from file without extension" << std::endl;
    return false;
  }

  if(ext == "obj")
    return read_OBJ(fname, points, polygons, np);
  else if(ext == "off")
    return read_OFF(fname, points, polygons, np);
  else if(ext == "ply")
    return read_PLY(fname, points, polygons, np);
  else if(ext == "stl")
    return read_STL(fname, points, polygons, np);
  else if(ext == "ts")
    return read_GOCAD(fname, points, polygons, np);
#ifdef CGAL_USE_VTK
  else if(ext == "ts")
    return read_VTP(fname, points, polygons, np);
#endif

  if(verbose)
  {
    std::cerr << "Error: unknown input file extension: " << ext << "\n"
              << "Please refer to the documentation for the list of supported file formats" << std::endl;
  }

  return false;
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange, typename PolygonRange>
bool read_polygon_soup(const std::string& fname, PointRange& points, PolygonRange& polygons)
{
  return read_polygon_soup(fname, points, polygons, parameters::all_default());
}

/// \endcond

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
 * \ingroup IOstreamFunctions
 *
 * \brief writes the content of `points` and `polygons` in a file.
 *
 * Supported file formats are the following:
 * - \ref IOStreamOFF (`.off`)
 * - \ref IOStreamOBJ (`.obj`)
 * - \ref IOStreamSTL (`.stl`)
 * - \ref IOStreamPLY (`.ply`)
 * - \ref IOStreamGocad (`.ts`)
 * - \ref IOStreamVTK (`.vtp`)
 *
 * The format is detected from the filename extension (letter case is not important).
 *
 * \tparam PolygonRange a model of the concept `RandomAccessContainer`
 *                      whose `value_type` is a model of the concept `RandomAccessContainer`
 *                      whose `value_type` is `std::size_t`.
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the name of the file.
 * \param points points of the soup of polygons
 * \param polygons each element in the range describes a polygon using the indices of the vertices.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{indicates whether output warnings and error messages should be printed or not.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if writing was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_polygon_soup(const std::string& fname,
                        const PointRange& points,
                        const PolygonRange& polygons,
                        const CGAL_BGL_NP_CLASS& np)
{
  const bool verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false);

  const std::string ext = internal::get_file_extension(fname);
  if(ext == std::string())
  {
    if(verbose)
      std::cerr << "Error: trying to output to file without extension" << std::endl;
    return false;
  }

  if(ext == "obj")
    return write_OBJ(fname, points, polygons, np);
  else if(ext == "off")
    return write_OFF(fname, points, polygons, np);
  else if(ext == "ply")
    return write_PLY(fname, points, polygons, np);
  else if(ext == "stl")
    return write_STL(fname, points, polygons, np);
  else if(ext == "ts")
    return write_GOCAD(fname, points, polygons, np);
#ifdef CGAL_USE_VTK
  else if(ext == "vtp")
    return write_VTP(fname, points, polygons, np);
#endif
#ifdef CGAL_LINKED_WITH_3MF
  else if(ext == "3mf")
    return write_3MF(fname, points, polygons);
#endif

  if(verbose)
  {
    std::cerr << "Error: unknown output file extension: " << ext << "\n"
              << "Please refer to the documentation for the list of supported file formats" << std::endl;
  }

  return false;
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange, typename PolygonRange>
bool write_polygon_soup(const std::string& fname, PointRange& points, PolygonRange& polygons)
{
  return write_polygon_soup(fname, points, polygons, parameters::all_default());
}

/// \endcond

} // namespace IO

} // namespace CGAL

#endif // CGAL_IO_READ_POLYGON_SOUP_H
