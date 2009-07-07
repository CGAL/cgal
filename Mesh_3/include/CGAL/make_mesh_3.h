// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description : make_mesh_3 function definition.
//******************************************************************************

#ifndef MAKE_MESH_3_H
#define MAKE_MESH_3_H

#include <CGAL/refine_mesh_3.h>

namespace CGAL {


/**
 * @brief This function meshes the domain defined by mesh_traits
 * (respecting criteria), and outputs the mesh to c3t3
 *
 * @param domain the domain to be discretized
 * @param criteria the criteria
 * @param exude if it is set to \c true, an exudation step will be done at
 *   the end of the Delaunay refinment process
 *
 * @return The mesh as a C3T3 object
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
C3T3 make_mesh_3(const MeshDomain&   domain,
                 const MeshCriteria& criteria,
                 bool exude = true)
{
  typedef typename MeshDomain::Point_3 Point_3;
  typedef typename MeshDomain::Index Index;
  typedef std::vector<std::pair<Point_3, Index> > Initial_points_vector;
  typedef typename Initial_points_vector::iterator Ipv_iterator;
  typedef typename C3T3::Vertex_handle Vertex_handle;
  
  // Mesh initialization : get some points and add them to the mesh
  Initial_points_vector initial_points;
  domain.construct_initial_points_object()(std::back_inserter(initial_points));
  C3T3 c3t3;
  
  // Insert points and set their index and dimension
  for ( Ipv_iterator it = initial_points.begin() ;
        it != initial_points.end() ;
        ++it )
  {
    Vertex_handle v = c3t3.triangulation().insert(it->first);
    c3t3.set_dimension(v,2); // by construction, points are on surface
    c3t3.set_index(v,it->second);
  }
  
  // Build mesher and launch refinement process
  // Don't reset c3t3 as we just created it
  refine_mesh_3(c3t3, domain, criteria, exude, false);
  
  return c3t3;
};


}  // end namespace CGAL


#endif // MAKE_MESH_3_H
