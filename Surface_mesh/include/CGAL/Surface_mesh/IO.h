//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
// Copyright (C) 2014 GeometryFactory
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_SURFACE_MESH_IO_H
#define CGAL_SURFACE_MESH_IO_H

#include <CGAL/license/Surface_mesh.h>

#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <CGAL/Surface_mesh/IO/3MF.h>
#include <CGAL/Surface_mesh/IO/OFF.h>
#include <CGAL/Surface_mesh/IO/PLY.h>

#include <CGAL/assertions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/use.h>

#include <boost/array.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <stdexcept>

namespace CGAL {

// @todo move that to read_polygon_mesh in BGL

/// \addtogroup PkgSurfaceMeshIO
///
/// I/O functionality for `Surface_mesh`. The top-level functions
/// `read_mesh()` and `write_mesh()` dispatch on the available readers
/// according to the file extension. Currently only `OFF` files are
/// supported.
///
/// @{

#if 0
/// Read a file into a `Surface_mesh`. The extension of the
/// filename determines which reader is used.
///
/// Mapping from extension to reader:
/// - off/OFF -> `read_OFF()`
///
/// @param mesh The mesh that should contain the input.
/// @param filename The name of the file to be read.
///
/// @return `true`, if reading succeeded, `false` otherwise
///
#endif
template <typename K>
bool read_mesh(Surface_mesh<K>& mesh, const std::string& filename)
{
  // clear mesh before reading from file
  mesh.clear();

  // extract file extension
  std::string::size_type dot(filename.rfind("."));
  if(dot == std::string::npos) return false;
  std::string ext = filename.substr(dot+1, filename.length()-dot-1);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  // extension determines reader
  if(ext == "off")
  {
    return read_OFF(filename, mesh);
  }

  // we didn't find a reader module
  return false;
}

#if 0
/// Write a `Surface_mesh` to a file. The extension of the
/// filename determines which writer is used.
///
/// Mapping from extension to writer:
/// - off/OFF -> `write_off()`
///
/// @param mesh The mesh to be written.
/// @param filename The name of the file to be written.
///
/// @return `true`, if writing succeeded, `false` otherwise
///
#endif
template <typename K>
bool write_mesh(const Surface_mesh<K>& mesh, const std::string& filename)
{
  // extract file extension
  std::string::size_type dot(filename.rfind("."));
  if(dot == std::string::npos) return false;
  std::string ext = filename.substr(dot+1, filename.length()-dot-1);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  // extension determines reader
  if(ext == "off")
  {
    return write_OFF(filename, mesh);
  }

  // we didn't find a writer module
  return false;
}

/// @}

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_IO_H
