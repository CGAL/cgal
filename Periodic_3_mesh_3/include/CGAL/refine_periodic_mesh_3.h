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
// $URL: svn+ssh://mbogdanov@scm.gforge.inria.fr/svn/cgal/trunk/Mesh_3/include/CGAL/refine_mesh_3.h $
// $Id: refine_mesh_3.h 60688 2011-01-10 15:43:22Z lrineau $
//
//
// Author(s)     : Mikhail Bogdanov
//
//******************************************************************************
// File Description : refine_mesh_3 function declaration and implementation.
//******************************************************************************

#ifndef CGAL_REFINE_PERIODIC_MESH_3_H
#define CGAL_REFINE_PERIODIC_MESH_3_H

#include <CGAL/refine_mesh_3.h>

#include <CGAL/Periodic_mesh_3/config.h>

namespace CGAL {

BOOST_PARAMETER_FUNCTION(
  (void),
  refine_periodic_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*) (criteria,*) ) // nondeduced
  (deduced 
    (optional
      (exude_param, (parameters::internal::Exude_options), parameters::no_exude()) // another default parameter distinct from Mesh_3
      (perturb_param, (parameters::internal::Perturb_options), parameters::no_perturb()) // another default parameter distinct from Mesh_3
      (odt_param, (parameters::internal::Odt_options), parameters::no_odt())
      (lloyd_param, (parameters::internal::Lloyd_options), parameters::no_lloyd())
      (reset_param, (parameters::Reset), parameters::reset_c3t3())
    )
  )
)
{
  // TODO: Can we call refine_mesh_3_impl?
  return refine_periodic_mesh_3_impl(c3t3,
                            domain,
                            criteria,
                            exude_param,
                            perturb_param,
                            odt_param,
                            lloyd_param,
                            reset_param() );
}
      
  
/**
 * @brief This function refines the mesh c3t3 wrt domain & criteria
 *
 * @param c3t3 the mesh to be refined.
 * @param domain the domain to be discretized
 * @param criteria the criteria
 * @param exude if \c true, an exudation step will be done at
 *   the end of the Delaunay refinement process
 * @param perturb if \c true, an explicit vertex perturbation step will be
 *   done at the end of refinement process
 * @param reset_c3t3 if \c true, a new C3T3 will be construct from param c3t3.
 *   The new c3t3 keeps only the vertices (as NON-weighted points with their
 *   dimension and Index) of the triangulation. That allows to refine a mesh
 *   which has been exuded.
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
void refine_periodic_mesh_3_impl(C3T3& c3t3,
                        const MeshDomain&   domain,
                        const MeshCriteria& criteria,
                        const parameters::internal::Exude_options& exude,
                        const parameters::internal::Perturb_options& perturb,
                        const parameters::internal::Odt_options& odt,
                        const parameters::internal::Lloyd_options& lloyd,
                        bool reset_c3t3)
{
  typedef Mesh_3::Mesher_3<C3T3, MeshCriteria, MeshDomain> Mesher;
  typedef typename C3T3::Triangulation::Geom_traits Gt;

  // Reset c3t3 (i.e. remove weights) if needed
  if ( reset_c3t3 )
  {
    C3T3 tmp_c3t3;
    std::for_each(c3t3.triangulation().finite_vertices_begin(),
                  c3t3.triangulation().finite_vertices_end(),
                  details::Insert_vertex_in_c3t3<C3T3>(tmp_c3t3));
    // TODO: corners and edges are not restored
    c3t3.swap(tmp_c3t3);
  }

  // Build mesher and launch refinement process
  Mesher mesher (c3t3, domain, criteria);
  double refine_time = mesher.refine_mesh();

  // try to project bad vertices
  //projection_of_external_points_of_surface(c3t3, domain);

  // Odt
  if ( odt )
  {
    std::cout << "odt" << std::endl;

    odt_optimize_mesh_3(c3t3,
                        domain,
                        parameters::time_limit = odt.time_limit(),
                        parameters::max_iteration_number = odt.max_iteration_number(),
                        parameters::convergence = odt.convergence(),
                        parameters::freeze_bound = odt.bound());
  }

  // Lloyd
  if ( lloyd )
  {
    std::cout << "lloyd" << std::endl;

    lloyd_optimize_mesh_3(c3t3,
                          domain,
                          parameters::time_limit = lloyd.time_limit(),
                          parameters::max_iteration_number = lloyd.max_iteration_number(),
                          parameters::convergence = lloyd.convergence(),
                          parameters::freeze_bound = lloyd.bound());
  }
  
  // Perturbation
  if ( perturb )
  {
    std::cout << "perturb" << std::endl;

    double perturb_time_limit = refine_time;

    if ( perturb.is_time_limit_set() )
      perturb_time_limit = perturb.time_limit();

    perturb_mesh_3(c3t3,
                   domain,
                   parameters::time_limit = perturb_time_limit,
                   parameters::sliver_bound = perturb.bound());
  }

#if 1 // regular stuff is not yet supported by periodic triangulations
  // Exudation
  if ( exude )
  {
    double exude_time_limit = refine_time;

    if ( exude.is_time_limit_set() )
      exude_time_limit = exude.time_limit();

    exude_mesh_3(c3t3,
                 parameters::time_limit = exude_time_limit,
                 parameters::sliver_bound = exude.bound());
  }
#endif //

}

}  // end namespace CGAL


#endif // CGAL_REFINE_PERIODIC_MESH_3_H
