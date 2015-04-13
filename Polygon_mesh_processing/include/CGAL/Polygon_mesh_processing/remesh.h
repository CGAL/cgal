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
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_REMESH_H
#define CGAL_POLYGON_MESH_PROCESSING_REMESH_H

#include <CGAL/Polygon_mesh_processing/internal/remesh_impl.h>

namespace CGAL {

namespace Polygon_mesh_processing {

/**
* \ingroup PkgPolygonMeshProcessing
* implements section 6.5.3 "Incremental remeshing" from the PMP book
*/
template<typename PolygonMesh>
void incremental_triangle_based_remeshing(PolygonMesh& pmesh,
                                          const double& target_edge_length,
                                          const unsigned int& nb_iterations = 10)
{
  double low  = 4./5. * target_edge_length;
  double high = 4./3. * target_edge_length;

  typename internal::Incremental_remesher<PolygonMesh> remesher(pmesh);
  for (unsigned int i = 0; i < nb_iterations; ++i)
  {
    remesher.split_long_edges(high);
    remesher.collapse_short_edges(low, high);
    remesher.equalize_valences();
    remesher.tangential_relaxation();
    remesher.project_to_surface();
  }
}

} //end namespace Polygon_mesh_processing
} //end namespace CGAL

#endif
