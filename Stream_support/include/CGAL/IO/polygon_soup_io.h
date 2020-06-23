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

// #include <CGAL/IO/3MF.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/IO/OFF.h>
// #include <CGAL/IO/OI.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/IO/STL.h>
// #include <CGAL/IO/VRML.h>
#include <CGAL/IO/VTK.h>
// #include <CGAL/IO/WKT.h>
#include <CGAL/IO/GOCAD.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {

/*!
 * \ingroup IOstreamFunctions
 *
 * \brief Reads a polygon soup from a file.
 *
 * Supported file formats are the following:
 * - \ref IOStreamOFF (`.off`)
 * - \ref IOStreamOBJ (`.obj`)
 * - \ref IOStreamSTL (`.stl`)
 * - \ref IOStreamPLY (`.ply`)
 * - \ref IOStreamGocad (`.ts`)
 * - \ref IOStreamVTK (`.vtp`)
 *
 * \tparam PolygonRange a model of the concept `RandomAccessContainer`
 * whose value_type is a model of the concept `RandomAccessContainer`
 * whose value_type is `std::size_t`.
 * \tparam PointRange a model of the concept `RandomAccessContainer`
 * whose value type is the point type
 *
 * \param fname the name of the file.
 * \param polygons each element in the range describes a polygon
 * using the indices of the vertices.
 * \param points points of the soup of polygons
 * \param verbose: if `true`, will output warnings and error messages.
 *
 * \return `true` if reading was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange>
bool read_polygon_soup(const std::string& fname,
                       PointRange& points,
                       PolygonRange& polygons,
                       const bool verbose = false)
{
  const std::string ext = IO::internal::get_file_extension(fname);
  if(ext == std::string())
  {
    if(verbose)
      std::cerr << "Error: cannot read from file without extension" << std::endl;
    return false;
  }

  if(ext == "obj")
    return read_OBJ(fname, points, polygons);
  else if(ext == "off")
    return read_OFF(fname, points, polygons);
  else if(ext == "ply")
    return read_PLY(fname, points, polygons);
  else if(ext == "stl")
    return read_STL(fname, points, polygons);
  else if(ext == "ts")
    return read_GOCAD(fname, points, polygons);
#ifdef CGAL_USE_VTK
  else if(ext == "ts")
    return read_VTP(fname, points, polygons);
#endif

  if(verbose)
  {
    std::cerr << "Error: unknown input file extension: " << ext << "\n"
              << "Please refer to the documentation for the list of supported file formats" << std::endl;
  }

  return false;
}

/*!
 * \ingroup IOstreamFunctions
 *
 * \brief Writes the content of `points` and `polygons` in a file.
 *
 * Supported file formats are the following:
 * - \ref IOStreamOFF (`.off`)
 * - \ref IOStreamOBJ (`.obj`)
 * - \ref IOStreamSTL (`.stl`)
 * - \ref IOStreamPLY (`.ply`)
 * - \ref IOStreamGocad (`.ts`)
 * - \ref IOStreamVTK (`.vtp`)
 *
 * \tparam PolygonRange a model of the concept `RandomAccessContainer`
 * whose value_type is a model of the concept `RandomAccessContainer`
 * whose value_type is `std::size_t`.
 * \tparam PointRange a model of the concept `RandomAccessContainer`
 * whose value type is the point type
 *
 * \param fname the name of the file.
 * \param polygons each element in the range describes a polygon
 * using the indices of the vertices.
 * \param points points of the soup of polygons
 * \param verbose: if `true`, will output warnings and error messages.
 *
 * \return `true` if writing was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange>
bool write_polygon_soup(const std::string& fname,
                        const PointRange& points,
                        const PolygonRange& polygons,
                        const bool verbose = false)
{
  const std::string ext = IO::internal::get_file_extension(fname);
  if(ext == std::string())
  {
    if(verbose)
      std::cerr << "Error: trying to output to file without extension" << std::endl;
    return false;
  }

  if(ext == "obj")
    return write_OBJ(fname, points, polygons);
  else if(ext == "off")
    return write_OFF(fname, points, polygons);
  else if(ext == "ply")
    return write_PLY(fname, points, polygons);
  else if(ext == "stl")
    return write_STL(fname, points, polygons);
  else if(ext == "ts")
    return write_GOCAD(fname, points, polygons);
#ifdef CGAL_USE_VTK
  else if(ext == "vtp")
    return write_VTP(fname, points, polygons);
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

} // namespace CGAL

#endif // CGAL_IO_READ_POLYGON_SOUP_H
