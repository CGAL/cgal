// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description : make_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_MAKE_MESH_3_H
#define CGAL_MAKE_MESH_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Mesh_3/config.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/tags.h>
#include <CGAL/Mesh_3/Protect_edges_sizing_field.h>
#include <CGAL/STL_Extension/internal/Has_features.h>
#include <CGAL/Mesh_3/C3T3_helpers.h>

#include <boost/mpl/has_xxx.hpp>
#include <type_traits>

#include <atomic>

namespace CGAL {

// -----------------------------------
// Initialize c3t3 stuff
// -----------------------------------
namespace Mesh_3 {
namespace internal {

template < typename C3T3, typename MeshDomain, typename InitialPointsGenerator >
void
add_points_from_generator(C3T3& c3t3,
                          const MeshDomain& domain,
                          const int nb_initial_points,
                          const InitialPointsGenerator& generator)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Vertex_handle Vertex_handle;
  typedef CGAL::Mesh_3::Triangulation_helpers<Tr> Th;

  struct PointDimIndex
  {
    typename Tr::Point m_wpt;
    int m_dim;
    typename MeshDomain::Index m_index;
  };

  const auto& cwp = c3t3.triangulation().geom_traits().construct_weighted_point_3_object();

  std::vector<PointDimIndex> initial_points;
  auto push_initial_point = [&](const auto& initial_pt)->void
    {
      using T = std::remove_cv_t < std::remove_reference_t<decltype(initial_pt)>>;
      if constexpr (std::tuple_size_v<T> == 3)
      {
        const auto& [weighted_pt, dim, index] = initial_pt;
        initial_points.push_back(PointDimIndex{ weighted_pt, dim, index });
      }
      else
      {
        const auto& [pt, index] = initial_pt;
        initial_points.push_back(PointDimIndex{ cwp(pt), 2, index });
      }
    };

  if (nb_initial_points > 0)
    generator(boost::make_function_output_iterator(push_initial_point), nb_initial_points);
  else
    generator(boost::make_function_output_iterator(push_initial_point));


  // Insert points and set their index and dimension
  for (const auto& [wpoint, dimension, index] : initial_points)
  {
    if(Th().inside_protecting_balls(c3t3.triangulation(), Vertex_handle(), wpoint.point()))
      continue;

    Vertex_handle v = c3t3.triangulation().insert(wpoint);

    // v could be null if point is hidden
    if ( v != Vertex_handle() )
    {
      c3t3.set_dimension(v,dimension);
      c3t3.set_index(v,index);
    }
  }
}

template<typename C3T3, typename MeshDomain>
bool
needs_more_init(C3T3& c3t3, const MeshDomain& domain)
{
  // If c3t3 initialization is not sufficient (may happen if
  // the user has not specified enough points ), add some surface points

  if (c3t3.triangulation().dimension() != 3)
    return true;
  else // dimension is 3 but it may not be enough
  {
    CGAL::Mesh_3::C3T3_helpers<C3T3, MeshDomain> helper(c3t3, domain);
    helper.update_restricted_facets();

    if (c3t3.number_of_facets() == 0) {
      return true;
    }
    else
    {
      helper.update_restricted_cells();
      if (c3t3.number_of_cells() == 0) {
        return true;
      }
    }
    return false;
  }
}

template < typename C3T3, typename MeshDomain, typename MeshCriteria, typename InitializationOptions>
void
init_c3t3(C3T3& c3t3, const MeshDomain& domain, const MeshCriteria&,
          const int nb_initial_points,
          const InitializationOptions& init_options)
{
  // 1st insert points from initial_points range and initial_points_generator
  add_points_from_generator(c3t3, domain, nb_initial_points, init_options);

  // If c3t3 initialization is not sufficient (may happen if
  // the user has not specified enough points ), add some surface points

  // use mesh domain's Construct_initial_points to complete initialization
  if(needs_more_init(c3t3, domain))
  {
    add_points_from_generator(c3t3, domain, nb_initial_points,
                              domain.construct_initial_points_object());
  }
}

template < typename EdgeCriteria >
struct Edge_criteria_sizing_field_wrapper
{
  typedef typename EdgeCriteria::Index    Index;
  typedef typename EdgeCriteria::FT       FT;
  typedef typename EdgeCriteria::Point_3  Point_3;

  Edge_criteria_sizing_field_wrapper(const EdgeCriteria& ec) : ec_(ec) {}
  FT operator()(const Point_3& p, const int dim, const Index& index) const
  { return ec_.sizing_field(p,dim,index); }

private:
  // No need to copy EdgeCriteria here
  const EdgeCriteria& ec_;
};

template < typename EdgeCriteria >
struct Edge_criteria_distance_field_wrapper
{
  typedef typename EdgeCriteria::Index    Index;
  typedef typename EdgeCriteria::FT       FT;
  typedef typename EdgeCriteria::Point_3  Point_3;

  Edge_criteria_distance_field_wrapper(const EdgeCriteria& ec) : ec_(ec) {}
  FT operator()(const Point_3& p, const int dim, const Index& index) const
  { return ec_.distance_field(p,dim,index); }

  private:
  // No need to copy EdgeCriteria here
  const EdgeCriteria& ec_;
};

template < typename C3T3, typename MeshDomain, typename MeshCriteria>
void init_c3t3_with_features(C3T3& c3t3,
                             const MeshDomain& domain,
                             const MeshCriteria& criteria,
                             bool nonlinear = false,
                             std::size_t maximal_number_of_vertices = 0,
                             Mesh_error_code* pointer_to_error_code = 0
#ifndef CGAL_NO_ATOMIC
                             , std::atomic<bool>* pointer_to_stop = 0
#endif
                             )
{
  typedef typename MeshCriteria::Edge_criteria Edge_criteria;
  typedef Edge_criteria_sizing_field_wrapper<Edge_criteria> Sizing_field;

  if (criteria.edge_criteria_object().has_distance_field())
  {
    typedef Edge_criteria_distance_field_wrapper<Edge_criteria> Distance_field;
    CGAL::Mesh_3::Protect_edges_sizing_field<C3T3,MeshDomain,Sizing_field,Distance_field>
      protect_edges(c3t3,
                    domain,
                    Sizing_field(criteria.edge_criteria_object()),
                    criteria.edge_criteria_object().min_length_bound(),
                    Distance_field(criteria.edge_criteria_object()),
                    maximal_number_of_vertices,
                    pointer_to_error_code
#ifndef CGAL_NO_ATOMIC
                    , pointer_to_stop
#endif
                    );
    protect_edges.set_nonlinear_growth_of_balls(nonlinear);

    protect_edges(true);
  }
  else
  {
    CGAL::Mesh_3::Protect_edges_sizing_field<C3T3,MeshDomain,Sizing_field>
      protect_edges(c3t3,
                    domain,
                    Sizing_field(criteria.edge_criteria_object()),
                    criteria.edge_criteria_object().min_length_bound(),
                    CGAL::Mesh_3::NoDistanceFunction(),
                    maximal_number_of_vertices,
                    pointer_to_error_code
#ifndef CGAL_NO_ATOMIC
                    , pointer_to_stop
#endif
                    );
    protect_edges.set_nonlinear_growth_of_balls(nonlinear);

    protect_edges(true);
  }
}

// This class is only used as base for specializations of C3t3_initializer
// when MeshDomain::Has_features is a valid type and is defined to CGAL::Tag_true
//
// Its purpose is to make the protection process virtual because Periodic_3_mesh_3
// handles sharp features differently and has its own 'init_c3t3_with_features()' function,
// but everything else is identical.
template < typename C3T3, typename MeshDomain, typename MeshCriteria>
struct C3t3_initializer_base
{
  virtual ~C3t3_initializer_base() { }

  // Not calling 'init_c3t3_with_features' directly to leave it as a free function
  // outside of the C3T3_initializer class
  virtual void
  initialize_features(C3T3& c3t3,
                      const MeshDomain& domain,
                      const MeshCriteria& criteria,
                      const parameters::internal::Mesh_3_options& mesh_options)
  {
    return Mesh_3::internal::init_c3t3_with_features
      (c3t3, domain, criteria,
       mesh_options.nonlinear_growth_of_balls,
       mesh_options.maximal_number_of_vertices,
       mesh_options.pointer_to_error_code
#ifndef CGAL_NO_ATOMIC
       , mesh_options.pointer_to_stop_atomic_boolean
#endif
       );
  }
};

// C3t3_initializer: initialize c3t3
template < typename C3T3,
           typename MeshDomain,
           typename MeshCriteria,
           bool MeshDomainHasHasFeatures,
           typename HasFeatures = int>
struct C3t3_initializer {};

// Partial specialization of C3t3_initializer
// Handle cases where MeshDomain::Has_features is not a valid type
template < typename C3T3, typename MD, typename MC, typename HasFeatures>
struct C3t3_initializer < C3T3, MD, MC, false, HasFeatures>
{
  typedef parameters::internal::Mesh_3_options Mesh_3_options;
  typedef parameters::internal::Initialization_options<MD, C3T3> Default_init_options;

  template<typename InitOptions = Default_init_options>
  void operator()(C3T3& c3t3,
                  const MD& domain,
                  const MC& criteria,
                  bool with_features,
                  Mesh_3_options mesh_options = Mesh_3_options(),
                  const InitOptions& init_options = Default_init_options())
  {
    if ( with_features )
    {
      std::cerr << "Warning: you requested a mesh with features from a domain"
                << " without features !" << std::endl;
    }

    init_c3t3(c3t3,domain,criteria,
              mesh_options.number_of_initial_points,
              init_options);
  }
};

// Partial specialization of C3t3_initializer
// Handles cases where MeshDomain::Has_features is a valid type
template < typename C3T3, typename MD, typename MC, typename HasFeatures>
struct C3t3_initializer < C3T3, MD, MC, true, HasFeatures>
{
  typedef parameters::internal::Mesh_3_options Mesh_3_options;
  typedef parameters::internal::Initialization_options<MD, C3T3> Default_init_options;

  template<typename InitOptions = Default_init_options>
  void operator()(C3T3& c3t3,
                  const MD& domain,
                  const MC& criteria,
                  bool with_features,
                  Mesh_3_options mesh_options = Mesh_3_options(),
                  const InitOptions& init_options = Default_init_options())
  {
    C3t3_initializer < C3T3, MD, MC, true, typename MD::Has_features >()
      (c3t3,domain,criteria,with_features,mesh_options,init_options);
  }
};

// Partial specialization of C3t3_initializer
// Handles cases where MeshDomain::Has_features is a valid type and is defined
// to CGAL::Tag_true
template < typename C3T3, typename MD, typename MC >
struct C3t3_initializer < C3T3, MD, MC, true, CGAL::Tag_true>
  : public C3t3_initializer_base < C3T3, MD, MC >
{
  virtual ~C3t3_initializer() { }

  typedef parameters::internal::Mesh_3_options Mesh_3_options;
  typedef parameters::internal::Initialization_options<MD, C3T3> Default_init_options;

  template<typename InitOptions = Default_init_options>
  void operator()(C3T3& c3t3,
                  const MD& domain,
                  const MC& criteria,
                  bool with_features,
                  Mesh_3_options mesh_options = Mesh_3_options(),
                  const InitOptions& init_options = Default_init_options())
  {
    if ( with_features ) {
      this->initialize_features(c3t3, domain, criteria,mesh_options);

      // If the initial points are not provided by the default generator,
      // there is no need to count the restricted facets and cells for now
      // because more vertices will be inserted anyway.
      // The check will be done later in init_c3t3()
      bool use_default_initializer = false;
      if constexpr (std::is_same_v<InitOptions, Default_init_options>) // check default type
      {
        use_default_initializer = init_options.is_default(); //check it also has no additional vertices
      }

      // If c3t3 initialization from features initialization
      // is not sufficient (may happen if there is only
      // a planar curve as feature for example), add some surface points.
      if (!use_default_initializer || CGAL::Mesh_3::internal::needs_more_init(c3t3, domain))
      {
        init_c3t3(c3t3, domain, criteria,
                  mesh_options.number_of_initial_points, init_options);
      }
    }
    else { init_c3t3(c3t3,domain,criteria,
                     mesh_options.number_of_initial_points, init_options); }
  }
};

// Partial specialization of C3t3_initializer
// Handles cases where MeshDomain::Has_features is a valid type and is defined
// to CGAL::Tag_false
template < typename C3T3, typename MD, typename MC >
struct C3t3_initializer < C3T3, MD, MC, true, CGAL::Tag_false>
{
  typedef parameters::internal::Mesh_3_options Mesh_3_options;
  typedef parameters::internal::Initialization_options<MD, C3T3> Default_init_options;

  template<typename InitOptions = Default_init_options>
  void operator()(C3T3& c3t3,
                  const MD& domain,
                  const MC& criteria,
                  bool with_features,
                  Mesh_3_options mesh_options = Mesh_3_options(),
                  const InitOptions& init_options = Default_init_options())
  {
    if ( with_features )
    {
      std::cerr << "Warning: you requested a mesh with features from a domain"
                << " without features !" << std::endl;
    }

    init_c3t3(c3t3,domain,criteria,
              mesh_options.number_of_initial_points, init_options);
  }
};

} // end namespace internal
} // end namespace Mesh_3


// -----------------------------------
// make_mesh_3 stuff
// -----------------------------------

/*!
 * \ingroup PkgMesh3Functions
 *
 * The function `make_mesh_3()` is a 3D
 * mesh generator. It produces simplicial meshes which discretize
 * 3D domains.
 *
 * The mesh generation algorithm is a Delaunay refinement process
 * followed by an optimization phase.
 * The criteria driving the Delaunay refinement
 * process may be tuned to achieve the user needs with respect to
 * the size of mesh elements, the accuracy of boundaries approximation,
 * etc.
 *
 * The optimization phase is a sequence of optimization processes,
 * amongst the following available optimizers: an ODT-smoothing,
 * a Lloyd-smoothing, a sliver perturber, and a sliver exuder.
 * Each optimization process
 * can be activated or not,
 * according to the user requirements
 * and available time.
 * By default, only the perturber and the exuder are activated.
 * Note that the benefits of the exuder will be lost if the mesh
 * is further refined afterward, and that ODT-smoothing, Lloyd-smoothing,
 * and sliver perturber should never be called after the sliver exuder.
 * In the case of further refinement, only the sliver exuder can be used.
 *
 * The function outputs the mesh to an object which provides iterators to
 * traverse the resulting mesh data structure or can be written to a file
 * (see \ref Mesh_3_section_examples ).
 *
 * \tparam C3T3 either a model of the concept `MeshComplex_3InTriangulation_3` or
 *              of `MeshComplexWithFeatures_3InTriangulation_3` if `MD`
 *              is a model of `MeshDomainWithFeatures_3`.
 *              The type `C3T3` is in particular required to provide a nested type
 *              `C3T3::Triangulation` for the 3D triangulation
 *              embedding the mesh. The vertex and cell base classes of the
 *              triangulation `C3T3::Triangulation` are required to be models of the
 *              concepts `MeshVertexBase_3` and `MeshCellBase_3` respectively.
 *
 * \tparam MD either a model of the concept `MeshDomain_3` or of
 *            `MeshDomainWithFeatures_3` if 0 and 1-dimensional features
 *            of the input complex have to be accurately represented in the mesh.
 *
 * \tparam MC either a model of the concept `MeshCriteria_3` or a model
 *            of `MeshCriteriaWithFeatures_3` if the domain has exposed features.
 *
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param domain the domain used to create the `c3t3` parameter. It is the sole link through which the domain
 *               to be discretized is known by the mesh generation algorithm.
 * \param criteria specifies the size and shape requirements for mesh tetrahedra
 *                 and surface facets. These criteria form the rules which drive
 *                 the refinement process. All mesh elements satisfy those criteria
 *                 at the end of the refinement process.
 *                 In addition, if the domain has features, the argument
 *                 `criteria` provides a sizing field to guide the discretization
 *                 of 1-dimensional exposed features.
 *
 *  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
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
 *                             <LI> `parameters::no_exude()`
 *                             <LI> `parameters::exude()`
 *                           </UL>}
 *     \cgalParamDefault{`parameters::exude()`}
 *   \cgalParamSectionEnd
 *   \cgalParamSectionBegin{Mesh initialization with a functor}
 *     \cgalParamDescription{an `InitialPointsGenerator` can optionally be provided to start the meshing process.
 *                           The following named parameter controls this option:
 *                           <UL>
 *                             <LI> `parameters::initial_points_generator()`
 *                           </UL>}
 *     \cgalParamDefault{the domain's `construct_initial_points_object()`
 *                       will be called for the points initialization.}
 *     \cgalParamExtra{If the generator does not generate enough points,
 *                    the domain's `construct_initial_points_object()` will be called.}
 *     \cgalParamExtra{If the parameter `parameters::initial_points()` is set,
 *                    the functor will be called after insertion of the points.}
 *   \cgalParamSectionEnd
 *   \cgalParamSectionBegin{Mesh initialization with points}
 *    \cgalParamDescription{a `Range` of initial points, represented as
 *                          tuple-like objects made of `tuple-like` objects of `<Weighted_point_3, int, Index>` can optionally
 *                          be provided to start the meshing process.
 *                          `Weighted_point_3` is the point's position and weight,
 *                          `int` is the dimension of the minimal dimension subcomplex on which
 *                          the point lies, and
 *                          `Index` is the corresponding subcomplex index.
 *                          The following named parameter controls this option:
 *                          <UL>
 *                            <LI> `parameters::initial_points()`
 *                          </UL>}
 *    \cgalParamDefault{`std::vector<std::tuple<Weighted_point_3, int, Index>>()`}
 *    \cgalParamExtra{If this parameter is set,
 *                    the domain's `construct_initial_points_object()` will be called
 *                    only if there is no facet in the restricted Delaunay triangulation
 *                    after points insertion.}
 *    \cgalParamExtra{If the parameter `parameters::initial_points_generator()` is set,
 *                    the points will be inserted before calling the functor.}
 *   \cgalParamSectionEnd
 * \cgalNamedParamsEnd
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
 * \sa `refine_mesh_3()`
 * \sa `exude_mesh_3()`
 * \sa `perturb_mesh_3()`
 * \sa `lloyd_optimize_mesh_3()`
 * \sa `odt_optimize_mesh_3()`
 */
template<typename C3T3, typename MeshDomain, typename MeshCriteria, typename CGAL_NP_TEMPLATE_PARAMETERS>
C3T3 make_mesh_3(const MeshDomain& domain, const MeshCriteria& criteria, const CGAL_NP_CLASS& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::get_parameter_reference;

    parameters::internal::Exude_options exude_param = choose_parameter(get_parameter(np, internal_np::exude_options_param), parameters::exude().v);
    parameters::internal::Perturb_options perturb_param = choose_parameter(get_parameter(np, internal_np::perturb_options_param), parameters::perturb().v);
    parameters::internal::Odt_options odt_param = choose_parameter(get_parameter(np, internal_np::odt_options_param), parameters::no_odt().v);
    parameters::internal::Lloyd_options lloyd_param = choose_parameter(get_parameter(np, internal_np::lloyd_options_param), parameters::no_lloyd().v);
    parameters::internal::Features_options features_param = choose_parameter(get_parameter(np, internal_np::features_options_param), parameters::features(domain).v);
    parameters::internal::Mesh_3_options mesh_options_param = choose_parameter(get_parameter(np, internal_np::mesh_param), parameters::internal::Mesh_3_options());
    parameters::internal::Manifold_options manifold_options_param = choose_parameter(get_parameter(np, internal_np::manifold_param), parameters::internal::Manifold_options());

    // range of initial points
    using Initial_point = std::pair<typename MeshDomain::Point_3, typename MeshDomain::Index>;
    using Initial_points_range_ref = typename internal_np::Lookup_named_param_def<internal_np::initial_points_param_t,
                                                                              CGAL_NP_CLASS,
                                                                              std::vector<Initial_point>>::reference;
    using Initial_points_range = std::remove_cv_t<std::remove_reference_t<Initial_points_range_ref>>;
    std::vector<Initial_point> empty_vec;
    Initial_points_range initial_points = choose_parameter(get_parameter_reference(np, internal_np::initial_points_param), empty_vec);

    // initial points generator
    using Initial_points_generator = typename internal_np::Lookup_named_param_def<internal_np::initial_points_generator_param_t,
                                                                                  CGAL_NP_CLASS,
                                                                                  typename MeshDomain::Construct_initial_points>::reference;
    auto default_generator = domain.construct_initial_points_object();
    Initial_points_generator initial_points_generator = choose_parameter(get_parameter(np, internal_np::initial_points_generator_param),
                                                                         default_generator);
    const parameters::internal::Initialization_options<MeshDomain, C3T3, Initial_points_generator, Initial_points_range>
      initial_points_gen_param(initial_points_generator, initial_points);

    C3T3 c3t3;
    make_mesh_3_impl(c3t3, domain, criteria,
            exude_param, perturb_param, odt_param, lloyd_param,
            features_param.features(), mesh_options_param,
            manifold_options_param,
            initial_points_gen_param);
    return c3t3;
}

#ifndef DOXYGEN_RUNNING
// Overload handling parameters passed with operator=
template<typename C3T3, typename MeshDomain, typename MeshCriteria,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
C3T3 make_mesh_3(const MeshDomain& domain, const MeshCriteria& criteria,
                 const CGAL_NP_CLASS_1&  np1,
                 const CGAL_NP_CLASS_2&  np2,
                 const NP& ... nps)
{
  return make_mesh_3<C3T3>(domain, criteria, internal_np::combine_named_parameters(np1, np2, nps...));
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
template<class C3T3, class MeshDomain, class MeshCriteria, class InitPtsGenerator, class InitPtsVec>
void make_mesh_3_impl(C3T3& c3t3,
                      const MeshDomain&   domain,
                      const MeshCriteria& criteria,
                      const parameters::internal::Exude_options& exude,
                      const parameters::internal::Perturb_options& perturb,
                      const parameters::internal::Odt_options& odt,
                      const parameters::internal::Lloyd_options& lloyd,
                      const bool with_features,
                      const parameters::internal::Mesh_3_options&
                        mesh_options = parameters::internal::Mesh_3_options(),
                      const parameters::internal::Manifold_options&
                        manifold_options = parameters::internal::Manifold_options(),
                      const parameters::internal::Initialization_options<MeshDomain, C3T3, InitPtsGenerator, InitPtsVec>&
                        initialization_options = parameters::internal::Initialization_options<MeshDomain, C3T3, InitPtsGenerator, InitPtsVec>())
{
#ifdef CGAL_MESH_3_INITIAL_POINTS_NO_RANDOM_SHOOTING
  CGAL::get_default_random() = CGAL::Random(0);
#endif

  // Initialize c3t3
  Mesh_3::internal::C3t3_initializer<
    C3T3,
    MeshDomain,
    MeshCriteria,
    ::CGAL::internal::has_Has_features<MeshDomain>::value,
    int>()(c3t3,
           domain,
           criteria,
           with_features,
           mesh_options,
           initialization_options);

  CGAL_assertion( c3t3.triangulation().dimension() >= 2 );

  // Build mesher and launch refinement process
  // Don't reset c3t3 as we just created it
  refine_mesh_3(c3t3, domain, criteria,
                parameters::exude_options=exude, parameters::perturb_options=perturb, parameters::odt_options=odt, parameters::lloyd_options= lloyd,
                parameters::no_reset_c3t3(), parameters::mesh_options= mesh_options,
                parameters::manifold_option= manifold_options);
}

#endif //DOXYGEN_RUNNING
} // end namespace CGAL

#endif // CGAL_MAKE_MESH_3_H
