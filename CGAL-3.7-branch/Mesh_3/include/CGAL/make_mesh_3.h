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

#ifndef CGAL_MAKE_MESH_3_H
#define CGAL_MAKE_MESH_3_H

#include <CGAL/Mesh_3/global_parameters.h>
#include <CGAL/refine_mesh_3.h>

namespace CGAL {

  
// Manual redirections
// boost::parameter can't handle make_mesh_3 return_type alone...
template <typename C3T3, typename MD, typename MC>
C3T3 make_mesh_3(const MD& md, const MC& mc)
{ 
  C3T3 c3t3;
  make_mesh_3_bp(c3t3,md,mc);
  return c3t3;
}
  
template <typename C3T3, typename MD, typename MC,
  typename Arg1>
C3T3 make_mesh_3(const MD& md, const MC& mc, const Arg1& a1)
{ 
  C3T3 c3t3;
  make_mesh_3_bp(c3t3,md,mc,a1);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC,
  typename Arg1, typename Arg2>
C3T3 make_mesh_3(const MD& md, const MC& mc, const Arg1& a1, const Arg2& a2)
{
  C3T3 c3t3; 
  make_mesh_3_bp(c3t3,md,mc,a1,a2);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC,
  typename Arg1, typename Arg2, typename Arg3>
C3T3 make_mesh_3(const MD& md, const MC& mc, const Arg1& a1, const Arg2& a2,
                 const Arg3& a3)
{ 
  C3T3 c3t3;
  make_mesh_3_bp(c3t3,md,mc,a1,a2,a3);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC,
  typename Arg1, typename Arg2, typename Arg3, typename Arg4>
C3T3 make_mesh_3(const MD& md, const MC& mc, const Arg1& a1, const Arg2& a2,
                 const Arg3& a3, const Arg4& a4)
{ 
  C3T3 c3t3;
  make_mesh_3_bp(c3t3,md,mc,a1,a2,a3,a4);
  return c3t3;
}
  
  
  

BOOST_PARAMETER_FUNCTION(
  (void),
  make_mesh_3_bp,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*) (criteria,*) ) // nondeduced
  (deduced 
    (optional
      (exude_param, (parameters::internal::Exude_options), parameters::exude())
      (perturb_param, (parameters::internal::Perturb_options), parameters::perturb())
      (odt_param, (parameters::internal::Odt_options), parameters::no_odt())
      (lloyd_param, (parameters::internal::Lloyd_options), parameters::no_lloyd())
    )
  )
)
{
  make_mesh_3_impl(c3t3, domain, criteria,
                   exude_param, perturb_param, odt_param, lloyd_param);
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
void make_mesh_3_impl(C3T3& c3t3,
                      const MeshDomain&   domain,
                      const MeshCriteria& criteria,
                      const parameters::internal::Exude_options& exude,
                      const parameters::internal::Perturb_options& perturb,
                      const parameters::internal::Odt_options& odt,
                      const parameters::internal::Lloyd_options& lloyd)
{
  typedef typename MeshDomain::Point_3 Point_3;
  typedef typename MeshDomain::Index Index;
  typedef std::vector<std::pair<Point_3, Index> > Initial_points_vector;
  typedef typename Initial_points_vector::iterator Ipv_iterator;
  typedef typename C3T3::Vertex_handle Vertex_handle;
  
  // Mesh initialization : get some points and add them to the mesh
  Initial_points_vector initial_points;
  domain.construct_initial_points_object()(std::back_inserter(initial_points));
  //C3T3 c3t3;
  
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
  refine_mesh_3(c3t3, domain, criteria,
                exude, perturb, odt, lloyd, parameters::no_reset_c3t3());
}


}  // end namespace CGAL


#endif // CGAL_MAKE_MESH_3_H
