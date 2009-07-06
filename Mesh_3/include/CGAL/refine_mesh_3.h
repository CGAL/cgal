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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : refine_mesh_3 function declaration and implementation.
//******************************************************************************

#ifndef REFINE_MESH_3_H_
#define REFINE_MESH_3_H_

#include <CGAL/Mesh_3/Slivers_exuder.h>
#include <CGAL/Mesh_3/Mesher_3.h>

namespace CGAL {

/**
 * @brief This function refines the mesh c3t3 wrt domain & criteria
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
void refine_mesh_3(C3T3& c3t3,
                   const MeshDomain&   domain,
                   const MeshCriteria& criteria,
                   bool exude = true)
{
  typedef Mesh_3::Mesher_3<C3T3, MeshCriteria, MeshDomain> Mesher;
  typedef typename C3T3::Triangulation::Geom_traits Gt;
  
  // Build mesher and launch refinement process
  Mesher mesher (c3t3, domain, criteria);
  mesher.refine_mesh();
  
  // Exudation
  if ( exude )
  {
    typedef Mesh_3::Min_dihedral_angle_criterion<Gt> Sliver_criterion;
    //typedef Mesh_3::Radius_radio_criterion<Gt> Sliver_criterion;
    typedef typename Mesh_3::Slivers_exuder<C3T3, Sliver_criterion> Exuder;
    
    Exuder exuder(c3t3);
  
#ifdef CGAL_MESH_3_VERBOSE
    exuder.print_stats();
#endif // CGAL_MESH_3_VERBOSE
  
    exuder.pump_vertices();
  
#ifdef CGAL_MESH_3_VERBOSE
    exuder.print_stats();
#endif // CGAL_MESH_3_VERBOSE
  }
  
};

}  // end namespace CGAL


#endif // REFINE_MESH_3_H_
