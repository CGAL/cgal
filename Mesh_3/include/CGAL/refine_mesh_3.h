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
// File Description :
//
//******************************************************************************

#ifndef REFINE_MESH_3_H_
#define REFINE_MESH_3_H_

#include <CGAL/Mesh_3/Mesher_3.h>

namespace CGAL {

/**
 * @brief This function refines the mesh c3t3 wrt domain & criteria
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
void refine_mesh_3(C3T3& c3t3,
                   const MeshDomain&   domain,
                   const MeshCriteria& criteria)
{
  typedef Mesh_3::Mesher_3<C3T3, MeshCriteria, MeshDomain>   Mesher;

  // Build mesher and launch refinement process
  Mesher mesher (c3t3, domain, criteria);
  mesher.refine_mesh();
};

}  // end namespace CGAL


#endif // REFINE_MESH_3_H_
