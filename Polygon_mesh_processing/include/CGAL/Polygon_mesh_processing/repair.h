// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot,
//                 Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_H
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/remove_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/remove_self_intersections.h>

#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {

/// \ingroup PMP_repairing_grp
/// removes the isolated vertices from any polygon mesh.
/// A vertex is considered isolated if it is not incident to any simplex
/// of higher dimension.
///
/// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
///
/// @param pmesh the polygon mesh to be repaired
///
/// @return number of removed isolated vertices
///
template <class PolygonMesh>
std::size_t remove_isolated_vertices(PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  std::vector<vertex_descriptor> to_be_removed;

  for(vertex_descriptor v : vertices(pmesh))
  {
    if(CGAL::halfedges_around_target(v, pmesh).first == CGAL::halfedges_around_target(v, pmesh).second)
      to_be_removed.push_back(v);
  }

  std::size_t nb_removed = to_be_removed.size();
  for(vertex_descriptor v : to_be_removed)
    remove_vertex(v, pmesh);

  return nb_removed;
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_H
