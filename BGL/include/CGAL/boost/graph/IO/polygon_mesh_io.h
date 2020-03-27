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

#ifndef CGAL_BGL_IO_POLYGON_MESH_IO_H
#define CGAL_BGL_IO_POLYGON_MESH_IO_H

#include <CGAL/boost/graph/IO/GOCAD.h>
#include <CGAL/boost/graph/IO/INP.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/boost/graph/IO/STL.h>
#include <CGAL/boost/graph/IO/VTK.h>
#include <CGAL/boost/graph/IO/WRL.h>

#include <algorithm>
#include <iostream>
#include <string>

namespace CGAL {

// @todo also need named parameters overload
template <typename FaceGraph>
bool read_polygon_mesh(const std::string& filename,
                       FaceGraph& g)
{
  std::string::size_type dot(filename.rfind("."));
  if(dot == std::string::npos)
    return false;

  std::string ext = filename.substr(dot+1, filename.length() - dot - 1);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  // extension determines reader
  if(ext == "obj")
    return read_OBJ(filename, g);
  else if(ext == "off") // @fixme coff, stoff, etc.
    return read_OFF(filename, g);
  else if(ext == "ply")
    return read_PLY(filename, g);
  else if(ext == "stl")
    return read_STL(filename, g);

  std::cerr << "Cannot open file with extension: " << ext << std::endl;

  return false;
}

template <typename FaceGraph>
bool write_polygon_mesh(const std::string& filename,
                        const FaceGraph& g)
{
  std::string::size_type dot(filename.rfind("."));
  if(dot == std::string::npos)
    return false;

  std::string ext = filename.substr(dot+1, filename.length()-dot-1);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  // extension determines writer
  if(ext == "obj")
    return write_OBJ(filename, g);
  else if(ext == "off") // @fixme coff, stoff, etc.
    return write_OFF(filename, g);
  else if(ext == "ply")
    return write_PLY(filename, g);
  else if(ext == "stl")
    return write_STL(filename, g);

  std::cerr << "Cannot open file with extension: " << ext << std::endl;

  return false;
}

} // namespace CGAL

#endif // CGAL_BGL_IO_POLYGON_MESH_IO_H
