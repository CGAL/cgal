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

#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>

#include <CGAL/Surface_mesh/IO/3MF.h>
#include <CGAL/Surface_mesh/IO/OFF.h>
#include <CGAL/Surface_mesh/IO/PLY.h>

#include <CGAL/boost/graph/io.h>

#include <string>

namespace CGAL {

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgSurfaceMeshIOFuncDeprecated
  \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::read_polygon_mesh()` should be used instead.
*/
template <typename K>
CGAL_DEPRECATED bool read_mesh(Surface_mesh<K>& sm, const std::string& filename)
{
  return IO::read_polygon_mesh(filename, sm);
}

/*!
  \ingroup PkgSurfaceMeshIOFuncDeprecated
  \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::write_polygon_mesh()` should be used instead.
*/
template <typename K>
CGAL_DEPRECATED bool write_mesh(const Surface_mesh<K>& mesh, const std::string& filename)
{
  return IO::write_polygon_mesh(filename, mesh);
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_IO_H
