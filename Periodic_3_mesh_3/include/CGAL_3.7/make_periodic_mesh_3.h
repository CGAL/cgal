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
// $URL: svn+ssh://mbogdanov@scm.gforge.inria.fr/svn/cgal/trunk/Mesh_3/include/CGAL/make_mesh_3.h $
// $Id: make_mesh_3.h 52705 2009-10-23 10:27:15Z stayeb $
//
//
// Author(s)     : Mikhail Bogdanov
//
//******************************************************************************
// File Description : make_periodic_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_MAKE_MESH_3_H
#define CGAL_MAKE_MESH_3_H

#include <CGAL/refine_periodic_mesh_3.h>

namespace CGAL {

template <typename C3T3, typename PMD, typename MC>
C3T3 make_periodic_mesh_3(const PMD& domain, const MC& mc)
{ 
  typedef typename C3T3::Triangulation Tr;

  C3T3 c3t3;
  Tr& tr = c3t3.triangulation();
  tr.set_domain(domain.periodic_cuboid());

  make_periodic_mesh_3_impl(c3t3, domain, mc);
  return c3t3;
}
  
/**
 * @brief This function meshes the domain defined by mesh_traits
 * (respecting criteria), and outputs the mesh to c3t3
 *
 * @param domain the domain to be discretized
 * @param criteria the criteria
 * @param exude if it is set to \c true, an exudation step will be done at
 *   the end of the Delaunay refinement process
 *
 * @return The mesh as a C3T3 object
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
void make_periodic_mesh_3_impl(C3T3& c3t3,
                      const MeshDomain&   domain,
                      const MeshCriteria& criteria)
{
  typedef typename MeshDomain::Point_3 Point_3;
  typedef typename MeshDomain::Index Index;
  typedef std::vector<std::pair<Point_3, Index> > Initial_points_vector;
  typedef typename Initial_points_vector::iterator Ipv_iterator;
  typedef typename C3T3::Vertex_handle Vertex_handle;
  
  // Mesh initialization : get some points and add them to the mesh
  Initial_points_vector initial_points;
  domain.construct_initial_points_object()(std::back_inserter(initial_points));
  
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
  refine_periodic_mesh_3(c3t3, domain, criteria, false);
};


}  // end namespace CGAL


#endif // CGAL_MAKE_MESH_3_H
