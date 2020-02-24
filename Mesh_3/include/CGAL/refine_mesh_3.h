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
#include <CGAL/boost/parameter.h>
#include <CGAL/Mesh_3/config.h>
#include <CGAL/Mesh_3/Dump_c3t3.h>
#include <CGAL/Mesh_3/Mesher_3.h>
#include <CGAL/Mesh_error_code.h>
#include <CGAL/optimize_mesh_3.h>
#include <CGAL/atomic.h>

#include <boost/parameter/preprocessor.hpp>

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
      typedef CGAL::cpp11::atomic<bool>* Pointer_to_stop_atomic_boolean_t;
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

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4003) // not enough actual parameters for macro
#endif

// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/boost/parameter.h>
CGAL_IGNORE_BOOST_PARAMETER_NAME_WARNINGS


// -----------------------------------
// Perturb
// -----------------------------------
BOOST_PARAMETER_FUNCTION((internal::Perturb_options), perturb, tag,
                         (optional (time_limit_, *, internal::undef_parameter )
                                   (sliver_bound_, *, default_values::perturb_sliver_bound )))
{
  internal::Perturb_options options(true);

  if ( internal::undef_parameter != time_limit_ )
    options.set_time_limit(time_limit_);

  options.set_bound(sliver_bound_);

  return options;
}

inline internal::Perturb_options no_perturb() { return internal::Perturb_options(false); }

// -----------------------------------
// Exude
// -----------------------------------
BOOST_PARAMETER_FUNCTION((internal::Exude_options), exude, tag,
                         (optional (time_limit_, *, internal::undef_parameter )
                                   (sliver_bound_, *, default_values::exude_sliver_bound )))
{
  internal::Exude_options options(true);

  if ( internal::undef_parameter != time_limit_ )
    options.set_time_limit(time_limit_);

  options.set_bound(sliver_bound_);

  return options;
}

inline internal::Exude_options no_exude() { return internal::Exude_options(false); }

// -----------------------------------
// Odt
// -----------------------------------
BOOST_PARAMETER_FUNCTION((internal::Odt_options), odt, tag,
                         (optional (time_limit_, *, 0 )
                                   (max_iteration_number_, *, 0 )
                                   (convergence_, *, default_values::odt_convergence_ratio )
                                   (freeze_bound_, *, default_values::odt_freeze_ratio )))
{
  internal::Odt_options options(true);

  options.set_time_limit(time_limit_);
  options.set_bound(freeze_bound_);
  options.set_convergence(convergence_);
  options.set_max_iteration_number(max_iteration_number_);

  return options;
}

inline internal::Odt_options no_odt() { return internal::Odt_options(false); }

// -----------------------------------
// Lloyd
// -----------------------------------
BOOST_PARAMETER_FUNCTION((internal::Lloyd_options), lloyd, tag,
                         (optional (time_limit_, *, 0 )
                                   (max_iteration_number_, *, 0 )
                                   (convergence_, *, default_values::lloyd_convergence_ratio )
                                   (freeze_bound_, *, default_values::lloyd_freeze_ratio )))
{
  internal::Lloyd_options options(true);

  options.set_time_limit(time_limit_);
  options.set_bound(freeze_bound_);
  options.set_convergence(convergence_);
  options.set_max_iteration_number(max_iteration_number_);

  return options;
}

inline internal::Lloyd_options no_lloyd() { return internal::Lloyd_options(false); }

// -----------------------------------
// Manifold options ------------------
// -----------------------------------
BOOST_PARAMETER_FUNCTION((internal::Manifold_options), manifold_options, tag,
                         (optional
                          (mesh_topology_, (int), -1)
                         )
                        )
{
  internal::Manifold_options options;
  options.mesh_topology = mesh_topology_;
  return options;
}

inline internal::Manifold_options manifold()
{
  return internal::Manifold_options(
          internal::Manifold_options::MANIFOLD);
}
inline internal::Manifold_options manifold_with_boundary()
{
  return internal::Manifold_options(
          internal::Manifold_options::MANIFOLD_WITH_BOUNDARY);
}
inline internal::Manifold_options non_manifold()
{
  return internal::Manifold_options(
          internal::Manifold_options::NON_MANIFOLD);
}

// -----------------------------------
// Mesh options
// -----------------------------------

// Undocumented Boost parameter for refine_mesh_3 and make_mesh_3.
// Allows to dump the mesh at given stage of the mesh generation
// algorithm.
BOOST_PARAMETER_FUNCTION((internal::Mesh_3_options), mesh_3_options, tag,
                         (optional
                          (dump_after_init_prefix_, (std::string), "" )
                          (dump_after_refine_surface_prefix_, (std::string), "" )
                          (dump_after_refine_prefix_, (std::string), "" )
                          (dump_after_glob_opt_prefix_, (std::string), "" )
                          (dump_after_perturb_prefix_, (std::string), "" )
                          (dump_after_exude_prefix_, (std::string), "" )
                          (number_of_initial_points_, (int), -1)
			  (maximal_number_of_vertices_, (std::size_t), 0)
                          (nonlinear_growth_of_balls_, (bool), false)
			  (pointer_to_error_code_, (Mesh_error_code*), ((Mesh_error_code*)0))
			  (pointer_to_stop_atomic_boolean_, (internal::Mesh_3_options::Pointer_to_stop_atomic_boolean_t), ((internal::Mesh_3_options::Pointer_to_stop_atomic_boolean_t)0))
                          )
                         )
{
  internal::Mesh_3_options options;

  options.dump_after_init_prefix=dump_after_init_prefix_;
  options.dump_after_refine_surface_prefix=dump_after_refine_surface_prefix_;
  options.dump_after_refine_prefix=dump_after_refine_prefix_;
  options.dump_after_glob_opt_prefix=dump_after_glob_opt_prefix_;
  options.dump_after_perturb_prefix=dump_after_perturb_prefix_;
  options.dump_after_exude_prefix=dump_after_exude_prefix_;
  options.number_of_initial_points=number_of_initial_points_;
  options.nonlinear_growth_of_balls = nonlinear_growth_of_balls_;
  options.maximal_number_of_vertices=maximal_number_of_vertices_;
  options.pointer_to_error_code=pointer_to_error_code_;
#ifndef CGAL_NO_ATOMIC
  options.pointer_to_stop_atomic_boolean=pointer_to_stop_atomic_boolean_;
#endif

  return options;
}

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

CGAL_PRAGMA_DIAG_POP

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

// -----------------------------------
// Reset_c3t3 (undocumented)
// -----------------------------------
  CGAL_BOOLEAN_PARAMETER(Reset,reset_c3t3,no_reset_c3t3)
  // CGAL_BOOLEAN_PARAMETER defined in <CGAL/boost/parameter.h>


// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/boost/parameter.h>
CGAL_IGNORE_BOOST_PARAMETER_NAME_WARNINGS

// -----------------------------------
// Parameters
// -----------------------------------
BOOST_PARAMETER_NAME( exude_param )
BOOST_PARAMETER_NAME( perturb_param )
BOOST_PARAMETER_NAME( odt_param )
BOOST_PARAMETER_NAME( lloyd_param )
BOOST_PARAMETER_NAME( reset_param )
BOOST_PARAMETER_NAME( mesh_options_param )
BOOST_PARAMETER_NAME( manifold_options_param )

CGAL_PRAGMA_DIAG_POP

} // end namespace parameters


  
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
  refine_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*) (criteria,*) ) // nondeduced
  (deduced
    (optional
      (exude_param, (parameters::internal::Exude_options), parameters::exude())
      (perturb_param, (parameters::internal::Perturb_options), parameters::perturb())
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

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_REFINE_MESH_3_H
