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

namespace IO {
namespace internal {

std::string get_file_extension(const std::string fname)
{
  std::string::size_type dot(fname.rfind("."));
  if(dot == std::string::npos)
    return std::string();

  std::string ext = fname.substr(dot+1, fname.length() - dot - 1);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  return ext;
}

} // namespace internal
} // namespace IO

/*!
 * \ingroup IOstreamFunctions
 *
 * \brief reads a polygon soup from a file.
 *
 * \tparam PolygonRange a model of the concept `RandomAccessContainer`
 * whose value_type is a model of the concept `RandomAccessContainer`
 * whose value_type is `std::size_t`.
 * \tparam PointRange a model of the concept `RandomAccessContainer`
 * whose value type is the point type
 *
 * \param fname the name of the file. Its extension must be one of the following :
 * `.off` (\ref IOStreamOFF "OFF file format") , `.obj` (\ref IOStreamOBJ "OBJ file format"),
 * `.stl` (\ref IOStreamSTL "STL file format"), `.ply` (\ref IOStreamPLY "PLY file format")
 * or `.ts`(\ref IOStreamGocad "GOCAD file format").
 * \param polygons each element in the range describes a polygon
 * using the indices of the vertices.
 * \param points points of the soup of polygons
 *
 * \return `true` if reading was successful, `false` otherwise.
 *
 * \see \ref IOStreamOFF
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
 * \brief writes a polygon soup in a file.
 *
 * \tparam PolygonRange a model of the concept `RandomAccessContainer`
 * whose value_type is a model of the concept `RandomAccessContainer`
 * whose value_type is `std::size_t`.
 * \tparam PointRange a model of the concept `RandomAccessContainer`
 * whose value type is the point type
 *
 * \param fname the name of the file. Its extension must be one of the following :
 * `.off` (\ref IOStreamOFF "OFF file format") , `.obj` (\ref IOStreamOBJ "OBJ file format"),
 * `.stl` (\ref IOStreamSTL "STL file format"), `.ply` (\ref IOStreamPLY "PLY file format")
 * or `.ts`(\ref IOStreamGocad "GOCAD file format").
 * \param polygons each element in the range describes a polygon
 * using the indices of the vertices.
 * \param points points of the soup of polygons
 *
 * \return `true` if writing was successful, `false` otherwise.
 *
 * \see \ref IOStreamOFF
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
  else if(ext == "ts")
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
