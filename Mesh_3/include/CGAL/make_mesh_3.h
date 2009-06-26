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
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
C3T3 make_mesh_3(const MeshDomain&   domain,
                 const MeshCriteria& criteria)
{
  typedef typename MeshDomain::Point_3 Point_3;
  typedef typename MeshDomain::Index Index;
  typedef std::vector<std::pair<Point_3, Index> > Initial_points_vector;

  // Mesh initialization : get some points and add them to the mesh
  Initial_points_vector initial_points;
  domain.construct_initial_points_object()(std::back_inserter(initial_points));
  C3T3 c3t3;
  c3t3.insert_surface_points(initial_points.begin(), initial_points.end());

  // Build mesher and launch refinement process
  refine_mesh_3(c3t3, domain, criteria);
  
  return c3t3;
};


}  // end namespace CGAL


#endif // MAKE_MESH_3_H
