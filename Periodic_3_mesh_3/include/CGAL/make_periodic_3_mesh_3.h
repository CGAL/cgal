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
// Author(s)     : Stéphane Tayeb,
//                 Mikhail Bogdanov,
//                 Mael Rouxel-Labbé
//
//******************************************************************************
// File Description : make_periodic_3_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_PERIODIC_3_MESH_3_MAKE_PERIODIC_3_MESH_3_H
#define CGAL_PERIODIC_3_MESH_3_MAKE_PERIODIC_3_MESH_3_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/Periodic_3_mesh_3/config.h>
#include <CGAL/Periodic_3_mesh_3/Protect_edges_sizing_field.h>
#include <CGAL/refine_periodic_3_mesh_3.h>

#include <CGAL/assertions.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_3/C3T3_helpers.h>
#include <CGAL/Named_function_parameters.h>

#include <boost/random/random_number_generator.hpp>
#include <boost/random/linear_congruential.hpp>

namespace CGAL {
namespace Periodic_3_mesh_3 {
namespace internal {

template<typename C3T3, typename MeshDomain>
bool insert_dummy_points(C3T3& c3t3,
                         const MeshDomain& domain)
{
  typedef typename C3T3::Triangulation::Vertex_handle                         Vertex_handle;

  typename C3T3::Triangulation::Geom_traits::Construct_point_3 cp =
    c3t3.triangulation().geom_traits().construct_point_3_object();

  std::vector<Vertex_handle> dummy_vertices = c3t3.triangulation().insert_generic_dummy_points();
  CGAL_postcondition(c3t3.triangulation().is_1_cover());

  // Abuse the multi cover function to check if cells have a small-enough orthoradius
  c3t3.triangulation().update_cover_data_after_converting_to_27_sheeted_covering();
  if(!c3t3.triangulation().can_be_converted_to_1_sheet())
  {
    std::cerr << "Error: dummy points do not create a 1-cover" << std::endl;
    CGAL_postcondition(false);
    return false;
  }

  for(Vertex_handle dvh : dummy_vertices)
  {
    dvh->info().is_dummy_vertex = true;

    // @fixme
    // The real in-dimension of a dummy vertex is impossible to evaluate: it can fall
    // on a surface or even a curve, but there is no way to know this because the concept (and models)
    // of MeshDomain_3 do not have a `Patch_index Do_intersect()(Point_3)`.
    //
    // It seems that setting the dimension wrongly, to a too-low dimension does not cause
    // problems, whereas on the other way is problematic: if setting it to '3' but the dummy vertex
    // is part of the restricted Delaunay, then the refinement will refine until the vertex
    // is no longer part of the restricted Delaunay, creating ugly spots in the mesh...
    c3t3.set_dimension(dvh, 0);

    // @fixme
    // This should be setting patch indices if the dummy vertex is not in dim 3.
    // Could maybe do it later, re-assigning indices of dummy vertices once the c3t3
    // has been scanned, but still before refinement starts.
    auto opt_si = domain.is_in_domain_object()(cp(c3t3.triangulation().point(dvh)));
    if(opt_si.has_value())
      c3t3.set_index(dvh, domain.index_from_subdomain_index(*opt_si));
    else
      c3t3.set_index(dvh, 0);
  }

  return true;
}

template <typename C3T3, typename MeshDomain, typename MeshCriteria>
void init_c3t3_with_features(C3T3& c3t3,
                             const MeshDomain& domain,
                             const MeshCriteria& criteria,
                             bool nonlinear = false)
{
  typedef typename MeshCriteria::Edge_criteria                                Edge_criteria;
  typedef Mesh_3::internal::Edge_criteria_sizing_field_wrapper<Edge_criteria> Sizing_field;

  CGAL::Periodic_3_mesh_3::Protect_edges_sizing_field<C3T3, MeshDomain, Sizing_field>
    protect_edges(c3t3, domain, Sizing_field(criteria.edge_criteria_object()));
  protect_edges.set_nonlinear_growth_of_balls(nonlinear);
  protect_edges(true);
}

// C3t3_initializer: initialize c3t3
template <typename C3T3,
          typename MeshDomain,
          typename MeshCriteria,
          bool MeshDomainHasHasFeatures,
          typename HasFeatures = int>
struct C3t3_initializer_base
  : public CGAL::Mesh_3::internal::C3t3_initializer<
      C3T3, MeshDomain, MeshCriteria, MeshDomainHasHasFeatures, HasFeatures>
{
  typedef CGAL::Mesh_3::internal::C3t3_initializer<
            C3T3, MeshDomain, MeshCriteria,
            MeshDomainHasHasFeatures, HasFeatures>              Base;

  void operator()(C3T3& c3t3,
                  const MeshDomain& domain,
                  const MeshCriteria& criteria,
                  bool with_features,
                  const parameters::internal::Mesh_3_options& mesh_options)
  {
    c3t3.triangulation().set_domain(domain.bounding_box());

    // Enforce 1-cover by adding dummy points
    if(!insert_dummy_points(c3t3, domain))
      return;

    // Call the basic initialization from c3t3, which handles features and
    // adds a bunch of points on the surface
    Base::operator()(c3t3, domain, criteria, with_features, mesh_options);
  }
};

template <typename C3T3,
          typename MeshDomain,
          typename MeshCriteria,
          bool MeshDomainHasHasFeatures,
          typename HasFeatures = int>
struct C3t3_initializer
  : public C3t3_initializer_base<C3T3, MeshDomain, MeshCriteria, MeshDomainHasHasFeatures, HasFeatures>
{
  typedef C3t3_initializer_base<C3T3, MeshDomain, MeshCriteria,
                                MeshDomainHasHasFeatures, HasFeatures> Base;

  void operator()(C3T3& c3t3, const MeshDomain& domain, const MeshCriteria& criteria,
                  bool with_features,
                  const parameters::internal::Mesh_3_options& mesh_options)
  {
    return Base::operator()(c3t3, domain, criteria, with_features, mesh_options);
  }
};

// Specialization when the mesh domain has 'Has_features'
template <typename C3T3,
          typename MeshDomain,
          typename MeshCriteria,
          typename HasFeatures>
struct C3t3_initializer<C3T3, MeshDomain, MeshCriteria, true, HasFeatures>
{
  void operator()(C3T3& c3t3, const MeshDomain& domain, const MeshCriteria& criteria,
                  bool with_features,
                  const parameters::internal::Mesh_3_options& mesh_options)
  {
    C3t3_initializer<C3T3, MeshDomain, MeshCriteria, true, typename MeshDomain::Has_features>()
        (c3t3, domain, criteria, with_features, mesh_options);
  }
};

// Specialization when the mesh domain has 'Has_features' and it's set to CGAL::Tag_true
template < typename C3T3,
           typename MeshDomain,
           typename MeshCriteria>
struct C3t3_initializer<C3T3, MeshDomain, MeshCriteria, true, CGAL::Tag_true>
  : public C3t3_initializer_base<C3T3, MeshDomain, MeshCriteria, true, CGAL::Tag_true>
{
  typedef C3t3_initializer_base<C3T3, MeshDomain, MeshCriteria, true, CGAL::Tag_true> Base;

  virtual ~C3t3_initializer() { }

  // this override will be used when initialize_features() is called, in make_mesh_3.h
  virtual void
  initialize_features(C3T3& c3t3,
                      const MeshDomain& domain,
                      const MeshCriteria& criteria,
                      const parameters::internal::Mesh_3_options& mesh_options)
  {
    return Periodic_3_mesh_3::internal::init_c3t3_with_features
      (c3t3, domain, criteria, mesh_options.nonlinear_growth_of_balls);
  }

  void operator()(C3T3& c3t3, const MeshDomain& domain, const MeshCriteria& criteria,
                  bool with_features,
                  const parameters::internal::Mesh_3_options& mesh_options)
  {
    return Base::operator()(c3t3, domain, criteria, with_features, mesh_options);
  }
};

} // namespace internal
} // namespace Periodic_3_mesh_3

// -----------------------------------
// make_periodic_3_mesh_3 stuff
// -----------------------------------

/*!
 * \ingroup PkgPeriodic3Mesh3Functions
 *
 * The function `make_periodic_3_mesh_3()` is a 3D periodic mesh generator.
 * It produces simplicial meshes which discretize 3D periodic domains.
 * The periodic mesh generation algorithm is a Delaunay refinement process
 * followed by an optimization phase. The criteria driving the Delaunay refinement
 * process may be tuned to achieve the user needs with respect to
 * the size of mesh elements, the accuracy of boundaries approximation,
 * etc.
 * The optimization phase is a sequence of optimization processes,
 * amongst the following available optimizers: an ODT-smoothing,
 * a Lloyd smoothing, a sliver perturber, and a sliver exuder.
 * Each optimization process can be activated or not, according to the user requirements
 * and available time.
 * By default, only the perturber and the exuder are activated.
 * Note that the benefits of the exuder will be lost if the mesh
 * is further refined afterward, and that ODT-smoothing, Lloyd-smoothing,
 * and sliver perturber should never be called after the sliver exuder.
 * In the case of further refinement, only the sliver exuder can be used.
 * The function outputs the mesh to an object which provides iterators to
 * traverse the resulting mesh data structure or can be written to a file
 * (see \ref Periodic_3_mesh_3_section_examples ).
 *
 *
 * \tparam C3T3 is required to be a model of
 * the concept `MeshComplex_3InTriangulation_3`. This is the return type.
 * The type `C3T3` is in particular required to provide a nested type
 * `C3T3::Triangulation` for the 3D triangulation
 * embedding the mesh. The vertex and cell base classes of the
 * triangulation `C3T3::Triangulation` are required to be models of the
 * concepts `MeshVertexBase_3` and `MeshCellBase_3`
 * respectively.
 *
 * \tparam MD is required to be a model of
 * the concept `MeshDomain_3`, or of the refined concept
 * `MeshDomainWithFeatures_3`
 * if the domain has corners and curves that need to be accurately represented in the mesh.
 * The argument `domain`
 * is the sole link through which the domain
 * to be discretized is known by the mesh generation algorithm.
 *
 * \tparam MC has to be a model of the concept
 * `MeshCriteria_3`, or a model of the refined concept `MeshCriteriaWithFeatures_3` if the domain has exposed features.
 * The argument `criteria` of type `MC` specifies the
 * size and shape requirements for mesh tetrahedra
 * and surface facets. These criteria
 * form the rules which drive the refinement process. All mesh elements
 * satisfy those criteria at the end of the refinement process.
 * In addition, if the domain has features, the argument
 * `criteria` provides a sizing field to guide the discretization
 * of 1-dimensional exposed features.
 *
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 *
 * \param domain the domain to be discretized
 * \param criteria the criteria
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamSectionBegin{Feature preservation options}
 *     \cgalParamDescription{If the domain is a model of `MeshDomainWithFeatures_3`, 0 and 1-dimensional features can be
 *                           taken into account while generating the mesh. The following two named parameters control
 *                           this option:
 *                           <UL>
 *                             <LI>\link parameters::features() `parameters::features(domain)` \endlink
 *                             <LI>`parameters::no_features()`
 *                           </UL>}
 *     \cgalParamDefault{`parameters::features(domain)`}
 *   \cgalParamSectionEnd
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
 *  The optimization parameters can be passed in an arbitrary order. If one parameter
 * is not passed, its default value is used. The default values are
 * `no_lloyd()`, `no_odt()`, `perturb()` and `exude()`.
 *
 * Note that regardless of which optimization processes are activated,
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
 * \sa `refine_periodic_3_mesh_3()`
 * \sa `make_mesh_3()`
 * \sa `exude_mesh_3()`
 * \sa `perturb_mesh_3()`
 * \sa `lloyd_optimize_mesh_3()`
 * \sa `odt_optimize_mesh_3()`
 */
template<typename C3T3, typename MeshDomain, typename MeshCriteria, typename CGAL_NP_TEMPLATE_PARAMETERS>
C3T3 make_periodic_3_mesh_3(const MeshDomain& domain, const MeshCriteria& criteria, const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  C3T3 c3t3;
  parameters::internal::Exude_options exude_param = choose_parameter(get_parameter(np, internal_np::exude_options_param), parameters::exude().v);
  parameters::internal::Perturb_options perturb_param = choose_parameter(get_parameter(np, internal_np::perturb_options_param), parameters::perturb().v);
  parameters::internal::Odt_options odt_param = choose_parameter(get_parameter(np, internal_np::odt_options_param), parameters::no_odt().v);
  parameters::internal::Lloyd_options lloyd_param = choose_parameter(get_parameter(np, internal_np::lloyd_options_param), parameters::no_lloyd().v);
  parameters::internal::Features_options features_param = choose_parameter(get_parameter(np, internal_np::features_options_param), parameters::features(domain).v);
  parameters::internal::Mesh_3_options mesh_options_param = choose_parameter(get_parameter(np, internal_np::mesh_param), parameters::internal::Mesh_3_options());
  parameters::internal::Manifold_options manifold_options_param = choose_parameter(get_parameter(np, internal_np::manifold_param), parameters::internal::Manifold_options());

  make_periodic_3_mesh_3_impl(c3t3, domain, criteria,
                   exude_param, perturb_param, odt_param, lloyd_param,
                   features_param.features(), mesh_options_param,
                   manifold_options_param);
  return c3t3;
}


#ifndef DOXYGEN_RUNNING
// Overload handling parameters passed with operator=
template<typename C3T3, typename MeshDomain, typename MeshCriteria,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
C3T3 make_periodic_3_mesh_3(const MeshDomain& domain, const MeshCriteria& criteria,
                            const CGAL_NP_CLASS_1&  np1,
                            const CGAL_NP_CLASS_2&  np2,
                            const NP& ... nps)
{
  return make_periodic_3_mesh_3<C3T3>(domain, criteria, internal_np::combine_named_parameters(np1, np2, nps...));
}


/**
 * @brief This function meshes the domain defined by mesh_traits
 * (respecting criteria), and outputs the mesh to c3t3
 *
 * @param domain the domain to be discretized
 * @param criteria the criteria
 * @param exude if it is set to `true`, an exudation step will be done at
 *   the end of the Delaunay refinement process
 *
 * @return The mesh as a C3T3 object
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
void make_periodic_3_mesh_3_impl(C3T3& c3t3,
                                 const MeshDomain& domain,
                                 const MeshCriteria& criteria,
                                 const parameters::internal::Exude_options& exude,
                                 const parameters::internal::Perturb_options& perturb,
                                 const parameters::internal::Odt_options& odt,
                                 const parameters::internal::Lloyd_options& lloyd,
                                 const bool with_features,
                                 const parameters::internal::Mesh_3_options&
                                   mesh_options = parameters::internal::Mesh_3_options(),
                                 const parameters::internal::Manifold_options&
                                   manifold_options = parameters::internal::Manifold_options())
{
  // Initialize c3t3
  Periodic_3_mesh_3::internal::C3t3_initializer<
    C3T3, MeshDomain, MeshCriteria,
    internal::has_Has_features<MeshDomain>::value>()(c3t3,
                                                     domain,
                                                     criteria,
                                                     with_features,
                                                     mesh_options);

  // Build mesher and launch refinement process
  refine_periodic_3_mesh_3(c3t3, domain, criteria,
                           parameters::exude_options = exude, parameters::perturb_options = perturb, parameters::odt_options = odt,
                           parameters::lloyd_options = lloyd, parameters::no_reset_c3t3(), // do not reset c3t3 as we just created it
                           parameters::mesh_options = mesh_options, parameters::manifold_option = manifold_options);
}
#endif //DOXYGEN_RUNNING
} // end namespace CGAL

#endif // CGAL_PERIODIC_3_MESH_3_MAKE_PERIODIC_3_MESH_3_H
