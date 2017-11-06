// Copyright (c) 2009, 2017 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mikhail Bogdanov
//
#ifndef CGAL_REFINE_PERIODIC_3_MESH_3_H
#define CGAL_REFINE_PERIODIC_3_MESH_3_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/Mesh_3/config.h>
#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/refine_mesh_3.h>
#include <CGAL/Mesh_3/Dump_c3t3.h>

#include <boost/parameter/preprocessor.hpp>

#include <algorithm>

namespace CGAL {

// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/Mesh_3/config.h>
CGAL_MESH_3_IGNORE_BOOST_PARAMETER_NAME_WARNINGS

BOOST_PARAMETER_FUNCTION(
  (void),
  refine_periodic_3_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*) (criteria,*) ) // nondeduced
  (deduced
    (optional
      (exude_param, (parameters::internal::Exude_options), parameters::no_exude()) // another default parameter distinct from Mesh_3
      (perturb_param, (parameters::internal::Perturb_options), parameters::no_perturb()) // another default parameter distinct from Mesh_3
      (odt_param, (parameters::internal::Odt_options), parameters::no_odt())
      (lloyd_param, (parameters::internal::Lloyd_options), parameters::no_lloyd())
      (reset_param, (parameters::Reset), parameters::reset_c3t3())
      (mesh_options_param, (parameters::internal::Mesh_3_options),
       parameters::internal::Mesh_3_options())
      (manifold_options_param, (parameters::internal::Manifold_options),
       parameters::internal::Manifold_options())
    )
  )
)
{
  return refine_periodic_3_mesh_3_impl(c3t3, domain, criteria,
                                       exude_param,
                                       perturb_param,
                                       odt_param,
                                       lloyd_param,
                                       reset_param(),
                                       mesh_options_param,
                                       manifold_options_param);
}

CGAL_PRAGMA_DIAG_POP

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
 * @param mesh_3_options is a struct object used to pass non-documented options,
 *   for debugging purpose.
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
void refine_periodic_3_mesh_3_impl(C3T3& c3t3,
                                   const MeshDomain& domain,
                                   const MeshCriteria& criteria,
                                   const parameters::internal::Exude_options& exude,
                                   const parameters::internal::Perturb_options& perturb,
                                   const parameters::internal::Odt_options& odt,
                                   const parameters::internal::Lloyd_options& lloyd,
                                   bool reset_c3t3,
                                   const parameters::internal::Mesh_3_options&
                                     mesh_options = parameters::internal::Mesh_3_options(),
                                   const parameters::internal::Manifold_options&
                                     manifold_options = parameters::internal::Manifold_options())
{
  typedef Mesh_3::Mesher_3<C3T3, MeshCriteria, MeshDomain> Mesher;

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

  dump_c3t3(c3t3, mesh_options.dump_after_init_prefix);

  // Build mesher and launch refinement process
  Mesher mesher(c3t3, domain, criteria,
                manifold_options.mesh_topology,
                mesh_options.maximal_number_of_vertices,
                mesh_options.pointer_to_error_code);

  double refine_time = mesher.refine_mesh(mesh_options.dump_after_refine_surface_prefix);
  c3t3.clear_manifold_info();

  // <periodic> below is the only difference with refine_mesh_3_impl. What does
  // it do and can it be removed ? (I guess something related to the initial grid)
  // try to project bad vertices
  //projection_of_external_points_of_surface(c3t3, domain);

  dump_c3t3(c3t3, mesh_options.dump_after_init_prefix);

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

  if( odt || lloyd) {
    dump_c3t3(c3t3, mesh_options.dump_after_glob_opt_prefix);
  }

  // Perturbation
  if ( perturb )
  {
    double perturb_time_limit = refine_time;

    if ( perturb.is_time_limit_set() )
      perturb_time_limit = perturb.time_limit();

    perturb_mesh_3(c3t3,
                   domain,
                   parameters::time_limit = perturb_time_limit,
                   parameters::sliver_bound = perturb.bound());

    dump_c3t3(c3t3, mesh_options.dump_after_perturb_prefix);
  }

  // Exudation
  if ( exude )
  {
    double exude_time_limit = refine_time;

    if ( exude.is_time_limit_set() )
      exude_time_limit = exude.time_limit();

    exude_mesh_3(c3t3,
                 parameters::time_limit = exude_time_limit,
                 parameters::sliver_bound = exude.bound());

    dump_c3t3(c3t3, mesh_options.dump_after_perturb_prefix);
  }
}

} // end namespace CGAL

#endif // CGAL_REFINE_PERIODIC_3_MESH_3_H
