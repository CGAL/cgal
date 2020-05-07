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
// #include <CGAL/IO/VTK.h>
// #include <CGAL/IO/WKT.h>
#include <CGAL/IO/GOCAD.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {

/*!
 * \ingroup IOstreamFunctions
 * \brief reads a polygon soup from a file.
 * \tparam PolygonRange a model of the concept `RandomAccessContainer`
 * whose value_type is a model of the concept `RandomAccessContainer`
 * whose value_type is `std::size_t`.
 * \tparam PointRange a model of the concept `RandomAccessContainer`
 * whose value type is the point type
 *
 * \param filename the name of the file. Its extension must be one of the following :
 * `.off` (\ref IOStreamOFF "OFF file format") , `.obj` (\ref IOStreamOBJ "OBJ file format"),
 * `.stl` (\ref IOStreamSTL "STL file format"), `.ply` (\ref IOStreamPLY "PLY file format")
 * or `.ts`(\ref IOStreamGocad "GOCAD file format").
 * \param polygons each element in the range describes a polygon
 * using the indices of the vertices.
 * \param points points of the soup of polygons

 * \return `true` if the reading worked, `false` otherwise.
 *
 * \see \ref IOStreamOFF
 */
template <typename PointRange, typename PolygonRange>
bool read_polygon_soup(const std::string& filename,
                       PointRange& points,
                       PolygonRange& polygons)
{
  std::string::size_type dot(filename.rfind("."));
  if(dot == std::string::npos)
    return false;

  std::string ext = filename.substr(dot+1, filename.length() - dot - 1);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  // extension determines reader
  if(ext == "obj")
    return read_OBJ(filename, points, polygons);
  else if(ext == "off")
    return read_OFF(filename, points, polygons);
  else if(ext == "ply")
    return read_PLY(filename, points, polygons);
  else if(ext == "stl")
    return read_STL(filename, points, polygons);
  else if(ext == "ts")
    return read_GOCAD(filename, points, polygons);

  std::cerr << "Cannot open file with extension: " << ext << std::endl;

  return false;
}

/*!
 * \ingroup IOstreamFunctions
 * \brief writes a polygon soup in a file.
 * \tparam PolygonRange a model of the concept `RandomAccessContainer`
 * whose value_type is a model of the concept `RandomAccessContainer`
 * whose value_type is `std::size_t`.
 * \tparam PointRange a model of the concept `RandomAccessContainer`
 * whose value type is the point type

 *
 * \param filename the name of the file. Its extension must be one of the following :
 * `.off` (\ref IOStreamOFF "OFF file format") , `.obj` (\ref IOStreamOBJ "OBJ file format"),
 * `.stl` (\ref IOStreamSTL "STL file format"), `.ply` (\ref IOStreamPLY "PLY file format")
 * or `.ts`(\ref IOStreamGocad "GOCAD file format").
 * \param polygons each element in the range describes a polygon
 * using the indices of the vertices.
 * \param points points of the soup of polygons

 * \return `true` if the writing worked, `false` otherwise.
 *
 * \see \ref IOStreamOFF
 */
template <typename PointRange, typename PolygonRange>
bool write_polygon_soup(const std::string& filename,
                        const PointRange& points,
                        const PolygonRange& polygons)
{
  std::string::size_type dot(filename.rfind("."));
  if(dot == std::string::npos)
    return false;

  std::string ext = filename.substr(dot+1, filename.length()-dot-1);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  // extension determines writer
  if(ext == "obj")
    return write_OBJ(filename, points, polygons);
  else if(ext == "off") // @fixme coff, stoff, etc.
    return write_OFF(filename, points, polygons);
  else if(ext == "ply")
    return write_PLY(filename, points, polygons);
  else if(ext == "stl")
    return write_STL(filename, points, polygons);
  else if(ext == "ts")
    return write_GOCAD(filename, points, polygons);

  std::cerr << "Cannot save file with extension: " << ext << std::endl;

  return false;
}

} // namespace CGAL

#endif // CGAL_IO_READ_POLYGON_SOUP_H
