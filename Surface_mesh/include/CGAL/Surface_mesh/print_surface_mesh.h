// Copyright (c) 2026 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Leo Valque
#ifndef CGAL_DEBUG_PRINT_SURFACE_MESH_H
#define CGAL_DEBUG_PRINT_SURFACE_MESH_H

#include <CGAL/license/Surface_mesh.h>

#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <iostream>

#ifndef DOXYGEN_RUNNING

namespace CGAL {

/// Print all mesh elements and their connectivity to an output stream.
template <typename P>
void print_surface_mesh(const Surface_mesh<P>& mesh,
                        std::ostream& out = std::cout)
{
  typedef Surface_mesh<P> Mesh;
  typedef typename Mesh::Vertex_index Vertex_index;
  typedef typename Mesh::Halfedge_index Halfedge_index;
  typedef typename Mesh::Face_index Face_index;

  out << "Points:\n";
  for (Vertex_index v : mesh.vertices())
    out << "  " << v << ": point=" << mesh.point(v) << '\n';

  out << "Vertices:\n";
  for (Vertex_index v : mesh.vertices())
    out << "  " << v << ": halfedge=" << mesh.halfedge(v) << '\n';

  out << "Halfedges:\n";
  for (Halfedge_index h : mesh.halfedges())
  {
    out << "  " << h << ": target=" << mesh.target(h)
        << " face=" << mesh.face(h)
        << " next=" << mesh.next(h)
        << " prev=" << mesh.prev(h) << '\n';
  }

  out << "Faces:\n";
  for (Face_index f : mesh.faces())
  {
    out << "  " << f << ": halfedge=" << mesh.halfedge(f) << '\n';
  }
}

} // namespace CGAL
#endif

#endif /* CGAL_DEBUG_PRINT_SURFACE_MESH_H */
