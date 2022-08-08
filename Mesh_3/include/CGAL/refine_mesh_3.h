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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : refine_mesh_3 function declaration and implementation.
//******************************************************************************

#ifndef CGAL_REFINE_MESH_3_H
#define CGAL_REFINE_MESH_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/config.h>
#include <CGAL/Mesh_3/config.h>
#include <CGAL/Mesh_3/Dump_c3t3.h>
#include <CGAL/Mesh_3/Mesher_3.h>
#include <CGAL/Mesh_error_code.h>
#include <CGAL/optimize_mesh_3.h>
#include <CGAL/Named_function_parameters.h>

#include <atomic>

#define CGAL_BOOLEAN_PARAMETER(Class, function_true, function_false)     \
  struct Class : public Base<bool> { Class(bool b) : Base<bool>(b){} };       \
  inline Named_function_parameters<Class, internal_np::reset_options_param_t> function_true() { \
  typedef Named_function_parameters<Class, internal_np::reset_options_param_t> Param;           \
  return Param(Class(true)); }                        \
  inline Named_function_parameters<Class, internal_np::reset_options_param_t> function_false() {\
  typedef Named_function_parameters<Class, internal_np::reset_options_param_t> Param;           \
  return Param(Class(false)); }


namespace CGAL {

namespace details {

/**
 * @class Insert_vertex_in_c3t3
 *
 * A functor designed to insert unweighted points into the triangulation
 * of a C3T3 from C3T3::Tr::Vertex , keeping the dimension and indices.
 */
template <typename C3T3>
class Insert_vertex_in_c3t3
{
private:
  typedef typename C3T3::Vertex_handle          Vertex_handle;
  typedef typename C3T3::Index                  Index;

  typedef typename C3T3::Triangulation          Tr;
  typedef typename Tr::Geom_traits              Geom_traits;
  typedef typename Tr::Vertex                   Vertex;
  typedef typename Tr::Weighted_point           Weighted_point;
  typedef typename Weighted_point::Weight       Weight;

public:
  Insert_vertex_in_c3t3(C3T3& c3t3)
    : r_c3t3_(c3t3) {}

  void operator()(const Vertex& vertex) const
  {
    typename Geom_traits::Construct_point_3 cp =
        r_c3t3_.triangulation().geom_traits().construct_point_3_object();
    typename Geom_traits::Compute_weight_3 cw =
        r_c3t3_.triangulation().geom_traits().compute_weight_3_object();

    // Get vh properties
    int dimension = vertex.in_dimension();
    Weight w = (dimension < 2) ? cw(vertex.point()) : 0;
    Weighted_point point(cp(vertex.point()), w);
    Index index = vertex.index();

    // Insert point and restore handle properties
    Vertex_handle new_vertex = r_c3t3_.triangulation().insert(point);
    r_c3t3_.set_index(new_vertex, index);
    r_c3t3_.set_dimension(new_vertex, dimension);

#if defined(CGAL_LINKED_WITH_TBB)\
&& !defined(CGAL_PARALLEL_MESH_3_DO_NOT_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE)
    if (boost::is_convertible<typename C3T3::Concurrency_tag, CGAL::Parallel_tag>::value)
    {
      if (dimension == -1)
        r_c3t3_.add_far_point(new_vertex);
    }
#endif
#ifdef CGAL_SEQUENTIAL_MESH_3_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE
    if (boost::is_convertible<typename C3T3::Concurrency_tag, CGAL::Sequential_tag>::value)
    {
      if (dimension == -1)
        r_c3t3_.add_far_point(new_vertex);
    }
#endif
  }

private:
  C3T3& r_c3t3_;
};

} // namespace details

namespace parameters {

namespace internal {

const int undef_parameter = -1;

// Helpers
struct Optimization_options_base
{
  Optimization_options_base(bool b)
  : b_(b), time_limit_(undef_parameter), bound_(undef_parameter) {}

  operator bool() const { return b_; }

  bool is_time_limit_set() const { return time_limit_ != undef_parameter; }
  void set_time_limit(double d) { time_limit_ = d; }
  double time_limit() const { return time_limit_; }

  bool is_bound_set() const { return bound_ != undef_parameter; }
  void set_bound(double d) { bound_ = d; }
  double bound() const { return bound_; }

private:
  bool b_;
  double time_limit_;
  double bound_;
};

struct Global_optimization_options_base
{
  Global_optimization_options_base()
  : convergence_(undef_parameter), max_it_nb_(undef_parameter) {}

  bool is_convergence_set() const { return convergence_ != undef_parameter; }
  void set_convergence(double d) { convergence_ = d; }
  double convergence() const { return convergence_; }

  bool is_max_iteration_number_set() const { return max_it_nb_ != undef_parameter; }
  void set_max_iteration_number(int i) { max_it_nb_ = i; }
  int max_iteration_number() const { return max_it_nb_; }

private:
  double convergence_;
  int max_it_nb_;
};

// Perturb
struct Perturb_options : public Optimization_options_base
{
  Perturb_options(bool b) : Optimization_options_base(b) {}
};

// Exude
struct Exude_options : public Optimization_options_base
{
  Exude_options(bool b) : Optimization_options_base(b) {}
};

// Odt
struct Odt_options : public Optimization_options_base
, public Global_optimization_options_base
{
  Odt_options(bool b) : Optimization_options_base(b)
  , Global_optimization_options_base() {}
};

// Lloyd
struct Lloyd_options : public Optimization_options_base
, public Global_optimization_options_base
{
  Lloyd_options(bool b) : Optimization_options_base(b)
  , Global_optimization_options_base() {}
};

// Manifold
struct Manifold_options {
  enum {
    NON_MANIFOLD = 0,
    MANIFOLD_WITH_BOUNDARY = 8,
    NO_BOUNDARY = 16,
    MANIFOLD = 24
  };

  Manifold_options(const int topology)
    : mesh_topology(topology)
  {}
  Manifold_options()
    : mesh_topology(NON_MANIFOLD)
  {}

  int mesh_topology;
};

// Various Mesh_3 option
struct Mesh_3_options {
#ifndef CGAL_NO_ATOMIC
      typedef std::atomic<bool>* Pointer_to_stop_atomic_boolean_t;
#else
      typedef bool* Pointer_to_stop_atomic_boolean_t;
#endif
  Mesh_3_options(bool nonlinear = false)
    // This parameter `nonlinear` adds a compatibility with previous
    // API of the constructor of `C3t3_initializer`.
    // -- Laurent Rineau, 2019/05/03
    : dump_after_init_prefix()
    , dump_after_refine_surface_prefix()
    , dump_after_refine_prefix()
    , dump_after_glob_opt_prefix()
    , dump_after_perturb_prefix()
    , dump_after_exude_prefix()
    , number_of_initial_points(-1)
    , nonlinear_growth_of_balls(nonlinear)
    , maximal_number_of_vertices(0)
    , pointer_to_error_code(0)
#ifndef CGAL_NO_ATOMIC
    , pointer_to_stop_atomic_boolean(0)
#endif
  {}

  std::string dump_after_init_prefix;
  std::string dump_after_refine_surface_prefix;
  std::string dump_after_refine_prefix;
  std::string dump_after_glob_opt_prefix;
  std::string dump_after_perturb_prefix;
  std::string dump_after_exude_prefix;
  int number_of_initial_points;
  bool nonlinear_growth_of_balls;
  std::size_t maximal_number_of_vertices;
  Mesh_error_code* pointer_to_error_code;
#ifndef CGAL_NO_ATOMIC
  Pointer_to_stop_atomic_boolean_t pointer_to_stop_atomic_boolean;
#endif

}; // end struct Mesh_3_options

} // end namespace internal

// Undocumented Boost parameter for refine_mesh_3 and make_mesh_3.
// Default Mesh_3_options: dump at every stage of the mesh generation.
inline internal::Mesh_3_options mesh_3_dump()
{
  internal::Mesh_3_options options;

  options.dump_after_init_prefix = "mesh_dump_after_init";
  options.dump_after_refine_surface_prefix = "mesh_dump_after_refine_surface";
  options.dump_after_refine_prefix = "mesh_dump_after_refine";
  options.dump_after_glob_opt_prefix = "mesh_dump_after_glob_opt";
  options.dump_after_perturb_prefix = "mesh_dump_after_perturb";
  options.dump_after_exude_prefix = "mesh_dump_after_exude";

  return options;
}
        template <typename T>
        struct Base
        {
            Base(T t) : t_(t) {}
            T operator()() const { return t_; }
        private:
            T t_;
        };

// -----------------------------------
// Reset_c3t3 (undocumented)
// -----------------------------------
  CGAL_BOOLEAN_PARAMETER(Reset,reset_c3t3,no_reset_c3t3)

} // end namespace parameters
#include <CGAL/Mesh_3/internal/parameters.h>
/*!
\ingroup PkgMesh3Functions

The function `refine_mesh_3()` is a 3D
mesh generator. It produces simplicial meshes which discretize
3D domains.

The mesh generation algorithm is a Delaunay refinement process
followed by an optimization phase.
The criteria driving the Delaunay refinement
process may be tuned to achieve the user needs with respect to
the size of mesh elements, the accuracy of boundaries approximation,
etc.

The optimization phase is a sequence of optimization processes,
amongst the following available optimizers: an ODT-smoothing,
a Lloyd smoothing, a sliver perturber, and a sliver exuder.
Each optimization process
can be activated or not,
according to the user requirements
and available time.
By default, only the perturber and the exuder are activated.
Note that the benefits of the exuder will be lost if the mesh
is further refined afterward.

\attention The function template `refine_mesh_3()` may be used to refine a previously
computed mesh, e.g.:
\code{.cpp}
C3T3 c3t3 = CGAL::make_mesh_3<C3T3>(domain,criteria);

CGAL::refine_mesh_3(c3t3, domain, new_criteria);
\endcode

Please note that we guarantee the result if and only if the domain does
not change from one refinement to the next one.


\tparam  C3T3 is required to be a model of
the concept
`MeshComplex_3InTriangulation_3`.
The argument `c3t3` is passed by
reference as this object is modified by the refinement process. As the
refinement process only adds points to the triangulation, all
vertices of the triangulation of `c3t3` remain in the
mesh during the refinement process. Object `c3t3` can be used to insert
specific points in the domain to ensure that they will be contained in the
final triangulation.
The type `C3T3` is in particular required to provide a nested type
`C3T3::Triangulation` for the 3D triangulation
embedding the mesh. The vertex and cell base classes of the
triangulation `C3T3::Triangulation` are required to be models of the
concepts `MeshVertexBase_3` and `MeshCellBase_3`
respectively.

\tparam MD is required to be a model of
the concept `MeshDomain_3` or of the refined concept
`MeshDomainWithFeatures_3` if 0 and 1-dimensional features
of the input complex have to be accurately represented in the mesh.
The argument `domain`
is the sole link through which the domain
to be discretized is known by the mesh generation algorithm.

\tparam MC is required to be a model of the concept
`MeshCriteria_3`, or a model of the refined concept `MeshCriteriaWithFeatures_3`
if the domain has exposed features. The argument `criteria` of
type `MC` specifies the
size and shape requirements for mesh tetrahedra
and surface facets. These criteria
form the rules which drive the refinement process. All mesh elements
satisfy those criteria at the end of the refinement process.
In addition, if the domain has features, the argument
`criteria` provides a sizing field to guide the discretization
of 1-dimensional exposed features.

\param c3t3 the mesh to be refined.
\param domain the domain to be discretized
\param criteria the criteria

The following four parameters are optional optimization parameters.
They control which optimization processes are performed
and allow the user to tune the parameters of the optimization processes.
Individual optimization parameters are not described here as they are
internal types (see instead the documentation page of each optimizer).
For each optimization algorithm, there exist two global functions
that allow to enable or disable the optimizer:

\cgalNamedParamsBegin
   \cgalParamNBegin{manifold_options_param}
     \cgalParamDescription{allows the user to drive the meshing algorithm,
                          and ensure that the output mesh surface follows the given manifold criterion.
                          It can be activated with `parameters::manifold()`, `parameters::manifold_with_boundary()`
                          and `parameters::non_manifold()`. Note that the meshing algorithm cannot generate a manifold
                          surface if the input surface is not manifold.}
      \cgalParamType{`parameters::manifold()` OR `parameters::manifold_with_boundary()` OR `parameters::non_manifold()`}
      \cgalParamDefault{`parameters::non_manifold()`}

    \cgalParamNBegin{lloyd_param}
     \cgalParamDescription{`parameters::lloyd()` and `parameters::no_lloyd()` are designed to
                            trigger or not a call to `lloyd_optimize_mesh_3()` function and to set the
                            parameters of this optimizer. If one parameter is not set, the default value of
                            `lloyd_optimize_mesh_3()` is used for this parameter.}
      \cgalParamType{`parameters::lloyd()` OR `parameters::no_lloyd()`}
      \cgalParamDefault{'parameters::no_lloyd()'}

    \cgalParamNBegin{odt_param}
     \cgalParamDescription{`parameters::odt()` and `parameters::no_odt()` are designed to
                            trigger or not a call to `odt_optimize_mesh_3()` function and
                            to set the parameters of this optimizer.
                            If one parameter is not set, the default value of
                           `odt_optimize_mesh_3()` is used for this parameter.}
      \cgalParamType{`parameters::odt()` OR `parameters::no_odt()`}
      \cgalParamDefault{`parameters::no_odt()`}

    \cgalParamNBegin{perturb_param}
     \cgalParamDescription{`parameters::perturb()` and `parameters::no_perturb()` are designed to
                            trigger or not a call to `perturb_mesh_3()` function and
                            to set the parameters of this optimizer. If one parameter is not set, the default value of
                            `perturb_mesh_3()` is used for this parameter, except for the time bound which is set to be
                             equal to the refinement CPU time.}
      \cgalParamType{`parameters::perturb()` and `parameters::no_perturb()`}
      \cgalParamDefault{'parameters::no_perturb`}

    \cgalParamNBegin{exude_param}
     \cgalParamDescription{parameters::exude()` and `parameters::no_exude()` are designed to
                           trigger or not a call to `exude_mesh_3()` function and to override to set the
                           parameters of this optimizer. If one parameter is not set, the default value of
                          `exude_mesh_3()` is used for this parameter, except for the time bound which is set to be
                           equal to the refinement CPU time.}
      \cgalParamType{`parameters::exude()` and `parameters::no_exude()`}
      \cgalParamDefault{'parameters::no_exude`}

\cgalNamedParamsEnd

The optimization parameters can be passed in arbitrary order. If one parameter
is not passed, its default value is used. The default values are
`no_lloyd()`, `no_odt()`, `perturb()` and `exude()`.
Note that whatever may be the optimization processes activated,
they are always launched in the order that is a suborder
of the following (see user manual for further
details): *ODT-smoother*, *Lloyd-smoother*, *perturber*, and *exuder*.

Beware that optimization of the mesh is obtained
by perturbing mesh vertices and modifying the mesh connectivity
and that this has an impact
on the strict compliance to the refinement criteria.
Though a strict compliance to mesh criteria
is guaranteed at the end of the Delaunay refinement, this may no longer be true after
some optimization processes. Also beware that the default behavior does involve some
optimization processes.

\sa `CGAL::make_mesh_3()`
\sa `CGAL::parameters::manifold`
\sa `CGAL::parameters::manifold_with_boundary`
\sa `CGAL::parameters::non_manifold`
\sa `CGAL::exude_mesh_3()`
\sa `CGAL::perturb_mesh_3()`
\sa `CGAL::lloyd_optimize_mesh_3()`
\sa `CGAL::odt_optimize_mesh_3()`
\sa `CGAL::parameters::exude`
\sa `CGAL::parameters::no_exude`
\sa `CGAL::parameters::perturb`
\sa `CGAL::parameters::no_perturb`
\sa `CGAL::parameters::lloyd`
\sa `CGAL::parameters::no_lloyd`
\sa `CGAL::parameters::odt`
\sa `CGAL::parameters::no_odt`
 */

template<typename C3T3, typename MeshDomain, typename MeshCriteria, typename CGAL_NP_TEMPLATE_PARAMETERS>
void refine_mesh_3(C3T3& c3t3, MeshDomain& domain, MeshCriteria& criteria, const CGAL_NP_CLASS& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    parameters::internal::Exude_options exude_param = choose_parameter(get_parameter(np, internal_np::exude_options_param), parameters::exude().v);
    parameters::internal::Perturb_options perturb_param = choose_parameter(get_parameter(np, internal_np::perturb_options_param), parameters::perturb().v);
    parameters::internal::Odt_options odt_param = choose_parameter(get_parameter(np, internal_np::odt_options_param), parameters::no_odt().v);
    parameters::internal::Lloyd_options lloyd_param = choose_parameter(get_parameter(np, internal_np::lloyd_options_param), parameters::no_lloyd().v);
    parameters::Reset reset_param = choose_parameter(get_parameter(np, internal_np::reset_options_param), parameters::reset_c3t3().v);
    parameters::internal::Mesh_3_options mesh_options_param = choose_parameter(get_parameter(np, internal_np::mesh_param), parameters::internal::Mesh_3_options());
    parameters::internal::Manifold_options manifold_options_param = choose_parameter(get_parameter(np, internal_np::manifold_param), parameters::internal::Manifold_options());

    return refine_mesh_3_impl(c3t3,
                              domain,
                              criteria,
                              exude_param,
                              perturb_param,
                              odt_param,
                              lloyd_param,
                              reset_param(),
                              mesh_options_param,
                              manifold_options_param);
}

#ifndef DOXYGEN_RUNNING
template<typename C3T3, typename MeshDomain, typename MeshCriteria, typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
void refine_mesh_3(C3T3& c3t3, MeshDomain& domain, MeshCriteria& criteria, const CGAL_NP_CLASS& ... nps)
{
    return refine_mesh_3(c3t3, domain, criteria, internal_np::combine_named_parameters(nps...));
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
 * @param mesh_3_options is a struct object used to pass non-documented options,
 *   for debugging purpose.
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
void refine_mesh_3_impl(C3T3& c3t3,
                        const MeshDomain&   domain,
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
  // Note: this function is almost entirely copied in refine_periodic_3_mesh.h
  // and any change to this function should likely be ported to the periodic version.

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
  Mesher mesher (c3t3, domain, criteria, manifold_options.mesh_topology,
                 mesh_options.maximal_number_of_vertices,
                 mesh_options.pointer_to_error_code
#ifndef CGAL_NO_ATOMIC
                 , mesh_options.pointer_to_stop_atomic_boolean
#endif
                 );
  double refine_time = mesher.refine_mesh(mesh_options.dump_after_refine_surface_prefix);
  c3t3.clear_manifold_info();

  dump_c3t3(c3t3, mesh_options.dump_after_refine_prefix);

  // Odt
  if ( odt )
  {
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

    dump_c3t3(c3t3, mesh_options.dump_after_exude_prefix);
  }
}
#endif // DOXYGEN_RUNNING
} // end namespace CGAL
#endif // CGAL_REFINE_MESH_3_H
