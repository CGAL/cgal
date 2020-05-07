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
