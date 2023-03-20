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

#include <CGAL/Periodic_3_mesh_3/config.h>
#include <CGAL/optimize_periodic_3_mesh_3.h>

#include <CGAL/Mesh_3/C3T3_helpers.h>
#include <CGAL/Mesh_3/Triangulation_helpers.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/SMDS_3/Dump_c3t3.h>
#include <CGAL/Time_stamper.h>

#include <CGAL/Named_function_parameters.h>
#include <boost/unordered_set.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>

namespace CGAL {
namespace internal {

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

    for(int i = 1; i < 4; i++)
    {
      Vertex_handle v = c->vertex((ind+i)&3);
      if(v->info().is_dummy_vertex)
      {
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_PROJECTION
        std::cout << c3t3.triangulation().point(v) << " must be projected" << std::endl;
#endif
        *vertices++ = v;
      }
    }
  }
}

template<class C3T3, class MeshDomain, class InputIterator>
bool project_points(C3T3& c3t3,
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
  typename C3T3::Triangulation::Geom_traits::Compute_squared_distance_3 csd =
    c3t3.triangulation().geom_traits().compute_squared_distance_3_object();

  CGAL::Mesh_3::C3T3_helpers<C3T3, MeshDomain> helper(c3t3, domain);
  CGAL::Mesh_3::Triangulation_helpers<Tr> tr_helpers;

  bool did_something = false;

  for(InputIterator it = vertex_begin; it != vertex_end; ++it)
  {
    Vertex_handle old_vertex = *it;

    const Weighted_point& weighted_old_position = c3t3.triangulation().point(old_vertex);
    CGAL_assertion(weighted_old_position.weight() == FT(0)); // point projection happens before optimizers

    const Bare_point& old_position = cp(weighted_old_position);
    const Bare_point new_position = helper.project_on_surface(old_vertex, old_position);
    const FT sq_d = csd(new_position, old_position);

#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_PROJECTION
    std::cerr << "\n\nMove dummy vertex" << std::endl;
    std::cerr << "old_vertex: " << &*old_vertex << std::endl;
    std::cerr << "old_position: " << old_position << std::endl;
    std::cerr << "new_position: " << new_position << std::endl;
    std::cerr << "squared distance from dummy to surface: " << sq_d << std::endl;
#endif

    // Skip tiny moves for efficiency
    auto min_v_and_sqd = c3t3.triangulation().nearest_power_vertex_with_sq_distance(old_vertex);
    CGAL_postcondition(min_v_and_sqd.first != Vertex_handle() && min_v_and_sqd.second != FT(-1));

    if(sq_d < 0.01 * min_v_and_sqd.second)
    {
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_PROJECTION
      std::cout << "REJECTED because dummy point is close enough to the surface" << std::endl;
#endif
      continue;
    }

    // Do not project if the projected point is in a protection ball
    if(tr_helpers.inside_protecting_balls(c3t3.triangulation(), old_vertex, new_position))
    {
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_PROJECTION
      std::cout << "REJECTED because new pos is within protection ball" << std::endl;
#endif
      continue;
    }

    // For periodic triangulations, the move is always performed using insert+remove,
    // so new_vertex cannot be old_vertex if the move has succeeded
    const Vector_3 move(old_position, new_position);
    Vertex_handle new_vertex = helper.update_mesh(old_vertex, move);

    // if the move has successfully been performed
    if(new_vertex != old_vertex && new_vertex != Vertex_handle())
    {
      new_vertex->info().is_dummy_vertex = false;
      c3t3.set_dimension(new_vertex, 2); // on the surface

      // @fixme
      // This actually should be the index from the surface patch index...
      // It can be obtained either by modifying project_on_surface to return the surface_patch index
      auto opt_si = domain.is_in_domain_object()(cp(c3t3.triangulation().point(new_vertex)));
      if(opt_si.has_value())
        c3t3.set_index(new_vertex, domain.index_from_subdomain_index(*opt_si));
      else
        c3t3.set_index(new_vertex, 0);
    }
    else
    {
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_PROJECTION
      std::cerr << "Warning: failed to create projection" << std::endl;
#endif
    }

    // The vertex `old_vertex` can still exist in the P3RT3:
    // - if the target already existed
    // - if its removal would have compromised the 1-cover property of the periodic triangulation
    // It's (almost) pointless to try and move it again, so fix it
    if(c3t3.triangulation().tds().is_vertex(old_vertex))
    {
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_PROJECTION
      std::cerr << "Warning: failed to remove pre-projection: " << c3t3.triangulation().point(old_vertex) << std::endl;
#endif
      old_vertex->info().is_dummy_vertex = false;
    }

    did_something = true;
  }

  return did_something;
}

template<class C3T3, class MeshDomain>
void project_dummy_points_of_surface(C3T3& c3t3,
                                     const MeshDomain& domain)
{
  typedef typename C3T3::Vertex_handle                     Vertex_handle;
  typedef CGAL::Hash_handles_with_or_without_timestamps    Hash_fct;
  typedef boost::unordered_set<Vertex_handle, Hash_fct>    Vertex_container;

  bool did_something = false;
  do
  {
    Vertex_container vertex_container;
    find_points_to_project(c3t3, std::insert_iterator<Vertex_container>(vertex_container, vertex_container.begin()));
    did_something = project_points(c3t3, domain, vertex_container.begin(), vertex_container.end());
  }
  while(did_something);
}

} // namespace internal

/*!
 * \ingroup PkgPeriodic3Mesh3Functions
 *
 * The function `refine_periodic_3_mesh_3()` is a 3D periodic
 * mesh generator. It produces periodic simplicial meshes which discretize
 * 3D periodic domains.
 *
 * The periodic mesh generation algorithm is a Delaunay refinement process
 * followed by an optimization phase.
 * The criteria driving the Delaunay refinement
 * process may be tuned to achieve the user needs with respect to
 * the size of mesh elements, the accuracy of boundaries approximation, etc.
 *
 * The optimization phase is a sequence of optimization processes,
 * amongst the following available optimizers: an ODT-smoothing,
 * a Lloyd smoothing, a sliver perturber, and a sliver exuder.
 * Each optimization process can be activated or not, according to the user requirements
 * and available time.
 * By default, only the perturber and the exuder are activated.
 * Note that the benefits of the exuder will be lost if the mesh
 * is further refined afterward.
 *
 * \attention The function template `refine_periodic_3_mesh_3()` may be used
 * to refine a previously computed mesh, e.g.:
 * \code{.cpp}
 * C3T3 c3t3 = CGAL::make_periodic_3_mesh_3<C3T3>(domain,criteria);
 *
 * CGAL::refine_periodic_3_mesh_3(c3t3, domain, new_criteria);
 * \endcode
 *
 * \attention Note that the triangulation must form at all times a simplicial complex within
 * a single copy of the domain (see Sections \ref P3Triangulation3secspace and \ref P3Triangulation3secintro
 * of the manual of 3D periodic triangulations). It is the responsibility of the user to provide
 * a triangulation that satisfies this condition when calling the refinement
 * function `refine_periodic_3_mesh_3`. The underlying triangulation of a mesh
 * complex obtained through `make_periodic_3_mesh_3()` or `refine_periodic_3_mesh_3()`
 * will always satisfy this condition.
 *
 *
 * \tparam  C3T3 is required to be a model of
 * the concept
 * `MeshComplex_3InTriangulation_3`.
 * The argument `c3t3` is passed by
 * reference as this object is modified by the refinement process. As the
 * refinement process only adds points to the triangulation, all
 * vertices of the triangulation of `c3t3` remain in the
 * mesh during the refinement process. Object `c3t3` can be used to insert
 * specific points in the domain to ensure that they will be contained in the
 * final triangulation.
 * The type `C3T3` is in particular required to provide a nested type
 * `C3T3::Triangulation` for the 3D triangulation
 * embedding the mesh. The vertex and cell base classes of the
 * triangulation `C3T3::Triangulation` are required to be models of the
 * concepts `MeshVertexBase_3` and `MeshCellBase_3`
 * respectively.
 *
 * \tparam MD is required to be a model of
 * the concept `Periodic_3MeshDomain_3` or of the refined concept
 * `Periodic_3MeshDomainWithFeatures_3` if 0 and 1-dimensional features
 * of the input complex have to be accurately represented in the mesh.
 * The argument `domain` is the sole link through which the domain
 * to be discretized is known by the mesh generation algorithm.
 *
 * \tparam MC is required to be a model of the concept
 * `MeshCriteria_3`, or a model of the refined concept `MeshCriteriaWithFeatures_3`
 * if the domain has exposed features. The argument `criteria` of
 * type `MC` specifies the
 * size and shape requirements for mesh tetrahedra
 * and surface facets. These criteria
 * form the rules which drive the refinement process. All mesh elements
 * satisfy those criteria at the end of the refinement process.
 * In addition, if the domain has features, the argument
 * `criteria` provides a sizing field to guide the discretization
 * of 1-dimensional exposed features.
 *
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters
 *
 * \param c3t3 the mesh to be refined.
 * \param domain the domain to be discretized
 * \param criteria the criteria
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 * The following four parameters are optional optimization parameters.
 * They control which optimization processes are performed
 * and allow the user to tune the parameters of the optimization processes.
 * Individual optimization parameters are not described here as they are
 * internal types (see instead the documentation page of each optimizer).
 * For each optimization algorithm, there exist two global functions
 * that allow to enable or disable the optimizer:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamSectionBegin{Topological options (manifoldness)}
 *     \cgalParamDescription{In order to drive the meshing algorithm and ensure that the output mesh follows a desired topological criterion,
 *                           three named parameters control this option:
 *                           <UL>
 *                             <LI>`parameters::manifold()`
 *                             <LI>`parameters::manifold_with_boundary()`
 *                             <LI>`parameters::non_manifold()`
 *                           </UL>
 *                           Note that the meshing algorithm cannot generate a manifold surface if the input surface is not manifold.}
 *     \cgalParamDefault{`parameters::non_manifold()`}
 *   \cgalParamSectionEnd
 *   \cgalParamSectionBegin{Lloyd optimization}
 *     \cgalParamDescription{`lloyd_optimize_mesh_3()` can optionally be called after the meshing process.
 *                           Two named parameters control this behavior:
 *                           <UL>
 *                             <LI> `parameters::no_lloyd()`
 *                             <LI> `parameters::lloyd_optimize_mesh_3()`
 *                           </UL>}
 *     \cgalParamDefault{`parameters::no_lloyd()`}
 *   \cgalParamSectionEnd
 *   \cgalParamSectionBegin{ODT optimization}
 *     \cgalParamDescription{`odt_optimize_mesh_3()` can optionally be called after the meshing process.
 *                           Two named parameters control this behavior:
 *                           <UL>
 *                             <LI> `parameters::no_odt()`
 *                             <LI> `parameters::odt()`
 *                           </UL>}
 *     \cgalParamDefault{`parameters::no_odt()`}
 *   \cgalParamSectionEnd
 *   \cgalParamSectionBegin{Mesh perturbation}
 *     \cgalParamDescription{`perturb_mesh_3()` can optionally be called after the meshing process.
 *                           Two named parameters control this behavior:
 *                           <UL>
 *                             <LI> `parameters::no_perturb()`
 *                             <LI> `parameters::perturb()`
 *                           </UL>}
 *     \cgalParamDefault{`parameters::perturb()`}
 *   \cgalParamSectionEnd
 *   \cgalParamSectionBegin{Mesh exudation}
 *     \cgalParamDescription{`exude_mesh_3()` can optionally be called after the meshing process.
 *                           Two named parameters control this behavior:
 *                           <UL>
 *                             <LI> `parameters::exude()`
 *                             <LI> `parameters::no_exude()`
 *                           </UL>}
 *     \cgalParamDefault{`parameters::exude()`}
 *   \cgalParamSectionEnd
 * \cgalNamedParamsEnd
 *
 * The optimization parameters can be passed in arbitrary order. If one parameter
 * is not passed, its default value is used. The default values are
 * `no_lloyd()`, `no_odt()`, `perturb()` and `exude()`.
 * Note that whatever may be the optimization processes activated,
 * they are always launched in the order that is a suborder
 * of the following (see user manual for further
 * details): *ODT-smoother*, *Lloyd-smoother*, *perturber*, and *exuder*.
 *
 * Beware that optimization of the mesh is obtained
 * by perturbing mesh vertices and modifying the mesh connectivity
 * and that this has an impact
 * on the strict compliance to the refinement criteria.
 * Though a strict compliance to mesh criteria
 * is guaranteed at the end of the Delaunay refinement, this may no longer be true after
 * some optimization processes. Also beware that the default behavior does involve some
 * optimization processes.
 *
 * \sa `make_periodic_3_mesh_3()`
 * \sa `refine_mesh_3()`
 * \sa `exude_periodic_3_mesh_3()`
 * \sa `perturb_periodic_3_mesh_3()`
 * \sa `lloyd_optimize_periodic_3_mesh_3()`
 * \sa `odt_optimize_periodic_3_mesh_3()`
 */
template<typename C3T3, typename MeshDomain, typename MeshCriteria, typename CGAL_NP_TEMPLATE_PARAMETERS>
void refine_periodic_3_mesh_3(C3T3& c3t3, const MeshDomain& domain, const MeshCriteria& criteria, const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  parameters::internal::Exude_options exude_param = choose_parameter(get_parameter(np, internal_np::exude_options_param), parameters::exude().v);
  parameters::internal::Perturb_options perturb_param = choose_parameter(get_parameter(np, internal_np::perturb_options_param), parameters::perturb().v);
  parameters::internal::Odt_options odt_param = choose_parameter(get_parameter(np, internal_np::odt_options_param), parameters::no_odt().v);
  parameters::internal::Lloyd_options lloyd_param = choose_parameter(get_parameter(np, internal_np::lloyd_options_param), parameters::no_lloyd().v);
  bool reset = choose_parameter(get_parameter(np, internal_np::do_reset_c3t3), false);
  parameters::internal::Mesh_3_options mesh_options_param = choose_parameter(get_parameter(np, internal_np::mesh_param), parameters::internal::Mesh_3_options());
  parameters::internal::Manifold_options manifold_options_param = choose_parameter(get_parameter(np, internal_np::manifold_param), parameters::internal::Manifold_options());

  return refine_periodic_3_mesh_3_impl(c3t3,
                                       domain,
                                       criteria,
                                       exude_param,
                                       perturb_param,
                                       odt_param,
                                       lloyd_param,
                                       reset,
                                       mesh_options_param,
                                       manifold_options_param);
}

#ifndef DOXYGEN_RUNNING
// Overload handling parameters passed with operator=
template<typename C3T3, typename MeshDomain, typename MeshCriteria,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
void refine_periodic_3_mesh_3(C3T3& c3t3, const MeshDomain& domain, const MeshCriteria& criteria,
                              const CGAL_NP_CLASS_1&  np1,
                              const CGAL_NP_CLASS_2&  np2,
                              const NP& ... nps)
{
  return refine_periodic_3_mesh_3(c3t3, domain, criteria, internal_np::combine_named_parameters(np1, np2, nps...));
}

/**
 * @brief This function refines the mesh c3t3 wrt domain & criteria
 *
 * @param c3t3 the mesh to be refined.
 * @param domain the domain to be discretized
 * @param criteria the criteria
 * @param exude if `true`, an exudation step will be done at
 *   the end of the Delaunay refinement process
 * @param perturb if `true`, an explicit vertex perturbation step will be
 *   done at the end of refinement process
 * @param reset_c3t3 if `true`, a new C3T3 will be construct from param c3t3.
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
    odt_optimize_periodic_3_mesh_3(c3t3, domain,
                                   parameters::time_limit = odt.time_limit(),
                                   parameters::max_iteration_number = odt.max_iteration_number(),
                                   parameters::convergence = odt.convergence(),
                                   parameters::freeze_bound = odt.bound());
  }

  // Lloyd
  if(lloyd)
  {
    lloyd_optimize_periodic_3_mesh_3(c3t3, domain,
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

    perturb_periodic_3_mesh_3(c3t3, domain,
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

    exude_periodic_3_mesh_3(c3t3,
                            parameters::time_limit = exude_time_limit,
                            parameters::sliver_bound = exude.bound());

    dump_c3t3(c3t3, mesh_options.dump_after_perturb_prefix);
  }

  CGAL_expensive_postcondition(c3t3.triangulation().tds().is_valid());
  CGAL_expensive_postcondition(c3t3.triangulation().is_valid());
  CGAL_expensive_postcondition(c3t3.is_valid());
}
#endif //DOXYGEN_RUNNING

} // end namespace CGAL

#endif // CGAL_REFINE_PERIODIC_3_MESH_3_H
