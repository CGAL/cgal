// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_PARAMETER_H
#define CGAL_BOOST_PARAMETER_H

#include <CGAL/config.h>


#ifdef BOOST_PARAMETER_MAX_ARITY
#  if (BOOST_PARAMETER_MAX_ARITY < 12)
#    error "BOOST_PARAMETER_MAX_ARITY must be at least 12 for CGAL::Mesh_3"
#  endif
#else
#  define  BOOST_PARAMETER_MAX_ARITY 12
#endif

#include <boost/parameter/name.hpp>

#if defined(__clang__) || (BOOST_GCC >= 40600)
#  define CGAL_IGNORE_UNUSED_VARIABLES \
    _Pragma("GCC diagnostic ignored \"-Wunused-variable\"") \
    _Pragma("GCC diagnostic ignored \"-Wunused-parameter\"")
#else
#  define CGAL_IGNORE_UNUSED_VARIABLES
#endif
#if __has_warning("-Wunneeded-internal-declaration")
#  define CGAL_IGNORE_UNUSED_INTERNAL_DECLARATION \
     _Pragma("clang diagnostic ignored \"-Wunneeded-internal-declaration\"")
#else
#  define CGAL_IGNORE_UNUSED_INTERNAL_DECLARATION
#endif

#define CGAL_IGNORE_BOOST_PARAMETER_NAME_WARNINGS \
  CGAL_IGNORE_UNUSED_VARIABLES                    \
  CGAL_IGNORE_UNUSED_INTERNAL_DECLARATION


namespace CGAL
{
namespace parameters
{

template <typename T>
struct Base
{
  Base(T t) : t_(t) {}
  T operator()() const { return t_; }
private:
  T t_;
};

#define CGAL_BOOLEAN_PARAMETER(Class, function_true, function_false)     \
  struct Class : public Base<bool> { Class(bool b) : Base<bool>(b){} };       \
  inline Class function_true() { return Class(true); }                        \
  inline Class function_false() { return Class(false); }

#define CGAL_DOUBLE_PARAMETER(Class, function, precondition)             \
  struct Class : public Base<double>                                          \
  { Class(double d) : Base<double>(d) { precondition(d); } };                 \
  inline Class function(double d) { return Class(d); }

// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/boost/parameter.h>
CGAL_IGNORE_BOOST_PARAMETER_NAME_WARNINGS

BOOST_PARAMETER_NAME( c3t3 )
BOOST_PARAMETER_NAME( domain )
BOOST_PARAMETER_NAME( criteria )
BOOST_PARAMETER_NAME( cdt )

BOOST_PARAMETER_NAME( (seeds_begin, tag) seeds_begin_)
BOOST_PARAMETER_NAME( (seeds_end, tag) seeds_end_)
BOOST_PARAMETER_NAME( (mark, tag) mark_)

BOOST_PARAMETER_NAME( (time_limit, tag) time_limit_ )
BOOST_PARAMETER_NAME( (convergence, tag) convergence_)
BOOST_PARAMETER_NAME( (max_iteration_number, tag) max_iteration_number_ )
BOOST_PARAMETER_NAME( (freeze_bound, tag) freeze_bound_)

BOOST_PARAMETER_NAME( (sliver_bound, tag) sliver_bound_)
BOOST_PARAMETER_NAME( (sliver_criterion, tag) sliver_criterion_)
BOOST_PARAMETER_NAME( (perturbation_vector, tag) perturbation_vector_)
BOOST_PARAMETER_NAME( (do_freeze, tag) do_freeze_)

BOOST_PARAMETER_NAME( (mesh_topology, tag) mesh_topology_)

BOOST_PARAMETER_NAME( (dump_after_init_prefix, tag ) dump_after_init_prefix_)
BOOST_PARAMETER_NAME( (dump_after_refine_surface_prefix, tag ) dump_after_refine_surface_prefix_)
BOOST_PARAMETER_NAME( (dump_after_refine_prefix, tag ) dump_after_refine_prefix_)
BOOST_PARAMETER_NAME( (dump_after_glob_opt_prefix, tag ) dump_after_glob_opt_prefix_)
BOOST_PARAMETER_NAME( (dump_after_perturb_prefix, tag ) dump_after_perturb_prefix_)
BOOST_PARAMETER_NAME( (dump_after_exude_prefix, tag ) dump_after_exude_prefix_)
BOOST_PARAMETER_NAME( (number_of_initial_points, tag) number_of_initial_points_)
BOOST_PARAMETER_NAME( (maximal_number_of_vertices, tag ) maximal_number_of_vertices_)
BOOST_PARAMETER_NAME( (nonlinear_growth_of_balls, tag ) nonlinear_growth_of_balls_)
BOOST_PARAMETER_NAME( (pointer_to_error_code, tag ) pointer_to_error_code_)
BOOST_PARAMETER_NAME( (pointer_to_stop_atomic_boolean, tag ) pointer_to_stop_atomic_boolean_)

// First used in <CGAL/Labeled_mesh_domain_3.h>
BOOST_PARAMETER_NAME( (function, tag ) function_)
BOOST_PARAMETER_NAME( (bounding_object, tag ) bounding_object_)
BOOST_PARAMETER_NAME( (relative_error_bound, tag ) relative_error_bound_)
BOOST_PARAMETER_NAME( (p_rng, tag ) p_rng_)
BOOST_PARAMETER_NAME( (null_subdomain_index, tag ) null_subdomain_index_)
BOOST_PARAMETER_NAME( (construct_surface_patch_index, tag ) construct_surface_patch_index_)

// First used in <CGAL/Gray_image_mesh_domain_3.h>
BOOST_PARAMETER_NAME( (image, tag ) image_)
BOOST_PARAMETER_NAME( (iso_value, tag) iso_value_)
BOOST_PARAMETER_NAME( (value_outside, tag) value_outside_)
BOOST_PARAMETER_NAME( (image_values_to_subdomain_indices, tag ) image_values_to_subdomain_indices_)

CGAL_PRAGMA_DIAG_POP
} // parameters
} // CGAL


#endif // CGAL_BOOST_PARAMETER_H
