// Copyright (c) 2009, 2014 INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Stephane Tayeb,
//                 Mikhail Bogdanov,
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_REFINE_PERIODIC_3_MESH_3_H
#define CGAL_REFINE_PERIODIC_3_MESH_3_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/Mesh_3/config.h>
#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Mesh_3/C3T3_helpers.h>
#include <CGAL/Mesh_3/Dump_c3t3.h>
#include <CGAL/Mesh_3/Triangulation_helpers.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Time_stamper.h>

#include <boost/parameter/preprocessor.hpp>
#include <boost/unordered_set.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>

namespace CGAL {

namespace internal {

template<class C3T3, class MeshDomain>
void project_dummy_points_of_surface(C3T3& c3t3, const MeshDomain& domain)
{
  typedef typename C3T3::Vertex_handle                     Vertex_handle;
  typedef CGAL::Hash_handles_with_or_without_timestamps    Hash_fct;
  typedef boost::unordered_set<Vertex_handle, Hash_fct>    Vertex_container;

  Vertex_container vertex_container;
  find_points_to_project(c3t3, std::insert_iterator<Vertex_container>(vertex_container, vertex_container.begin()));

  project_points(c3t3, domain, vertex_container.begin(), vertex_container.end());
}

template<class C3T3, class OutputIterator>
void find_points_to_project(C3T3& c3t3, OutputIterator vertices)
{
  typedef typename C3T3::Vertex_handle Vertex_handle;
  typedef typename C3T3::Cell_handle Cell_handle;

  // Don't project all dummy points, but simply the ones that are involved in
  // surface facets of the c3t3.
  for(typename C3T3::Facets_in_complex_iterator face_it = c3t3.facets_in_complex_begin();
                                                face_it != c3t3.facets_in_complex_end();
                                                ++face_it)
  {
    int ind = face_it->second;
    Cell_handle c = face_it->first;

    for(int i = 1; i < 4; i++) {
      Vertex_handle v = c->vertex((ind+i)&3);

      typename C3T3::Index index = c3t3.index(v);
      if(const int* i = boost::get<int>(&index))
      {
        if(*i == 0) // '0' is the index of dummies
          *vertices++ = v;
      }
    }
  }
}

template<class C3T3, class MeshDomain, class InputIterator>
void project_points(C3T3& c3t3,
                    const MeshDomain& domain,
                    InputIterator vertex_begin,
                    InputIterator vertex_end)
{
  typedef typename C3T3::Vertex_handle         Vertex_handle;

  typedef typename C3T3::Triangulation         Tr;

  typedef typename Tr::Bare_point              Bare_point;
  typedef typename Tr::Weighted_point          Weighted_point;
  typedef typename Tr::Geom_traits::Vector_3   Vector_3;
  typedef typename Tr::Geom_traits::FT         FT;

  typename C3T3::Triangulation::Geom_traits::Construct_point_3 cp =
    c3t3.triangulation().geom_traits().construct_point_3_object();

  CGAL::Mesh_3::C3T3_helpers<C3T3, MeshDomain> helper(c3t3, domain);
  CGAL::Mesh_3::Triangulation_helpers<Tr> tr_helpers;

  for(InputIterator it = vertex_begin; it != vertex_end; ++it)
  {
    Vertex_handle vh = *it;

    const Weighted_point& vh_wp = c3t3.triangulation().point(vh);
    const Bare_point& vh_p = cp(vh_wp);
    const Bare_point new_point = helper.project_on_surface(vh, vh_p);

    const FT sq_d = CGAL::squared_distance(new_point, vh_p);

#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_PROJECTION
    std::cerr << "vh: " << &*vh << std::endl;
    std::cerr << "vhp: " << vh_p << std::endl;
    std::cerr << "projected: " << new_point << std::endl;
    std::cerr << "squared distance from dummy to surface: " << sq_d << std::endl;
#endif

    // Skip tiny moves for efficiency
    if(sq_d < 1e-10) // arbitrary value, maybe compare it to the surface distance criterium ?
      continue;

    // Do not project if the projected point is in a protection ball
    if(tr_helpers.inside_protecting_balls(c3t3.triangulation(), vh, new_point))
      continue;

    const Vector_3 move(vh_p, new_point);
    Vertex_handle new_vertex = helper.update_mesh(vh, move);
    if(new_vertex != vh) // if the move has successfully been performed
      c3t3.set_dimension(new_vertex, 2);
  }
}

} // namespace internal

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4003) // not enough actual parameters for macro
#endif

// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/boost/parameter.h>
CGAL_IGNORE_BOOST_PARAMETER_NAME_WARNINGS

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

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

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
  // 'refine_periodic_3_mesh_3' does not insert dummy points, so the triangulation
  // must already be 1-sheeted prior to calls to this function.
  CGAL_precondition(c3t3.triangulation().is_1_cover());

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

  CGAL_expensive_postcondition(c3t3.triangulation().tds().is_valid());
  CGAL_expensive_postcondition(c3t3.triangulation().is_valid());
  CGAL_expensive_postcondition(c3t3.is_valid());
  dump_c3t3(c3t3, mesh_options.dump_after_init_prefix);

  // Project dummy points on the surface to lessen their influence on the output
  internal::project_dummy_points_of_surface(c3t3, domain);

  CGAL_expensive_postcondition(c3t3.triangulation().tds().is_valid());
  CGAL_expensive_postcondition(c3t3.triangulation().is_valid());
  CGAL_expensive_postcondition(c3t3.is_valid());

  // Odt
  if(odt)
  {
    odt_optimize_mesh_3(c3t3, domain,
                        parameters::time_limit = odt.time_limit(),
                        parameters::max_iteration_number = odt.max_iteration_number(),
                        parameters::convergence = odt.convergence(),
                        parameters::freeze_bound = odt.bound());
  }

  // Lloyd
  if(lloyd)
  {
    lloyd_optimize_mesh_3(c3t3, domain,
                          parameters::time_limit = lloyd.time_limit(),
                          parameters::max_iteration_number = lloyd.max_iteration_number(),
                          parameters::convergence = lloyd.convergence(),
                          parameters::freeze_bound = lloyd.bound());
  }

  if(odt || lloyd)
  {
    dump_c3t3(c3t3, mesh_options.dump_after_glob_opt_prefix);
  }

  // Perturbation
  if(perturb)
  {
    double perturb_time_limit = refine_time;

    if(perturb.is_time_limit_set())
      perturb_time_limit = perturb.time_limit();

    perturb_mesh_3(c3t3, domain,
                   parameters::time_limit = perturb_time_limit,
                   parameters::sliver_bound = perturb.bound());

    dump_c3t3(c3t3, mesh_options.dump_after_perturb_prefix);
  }

  // Exudation
  if(exude)
  {
    double exude_time_limit = refine_time;

    if(exude.is_time_limit_set())
      exude_time_limit = exude.time_limit();

    exude_mesh_3(c3t3,
                 parameters::time_limit = exude_time_limit,
                 parameters::sliver_bound = exude.bound());

    dump_c3t3(c3t3, mesh_options.dump_after_perturb_prefix);
  }

  CGAL_expensive_postcondition(c3t3.triangulation().tds().is_valid());
  CGAL_expensive_postcondition(c3t3.triangulation().is_valid());
  CGAL_expensive_postcondition(c3t3.is_valid());
}

} // end namespace CGAL

#endif // CGAL_REFINE_PERIODIC_3_MESH_3_H
