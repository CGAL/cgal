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
//
//
//******************************************************************************
// File Description : Named Function Parameters specific code
//******************************************************************************
#ifndef CGAL_PARAMETERS_H
#define CGAL_PARAMETERS_H

#include <CGAL/license/Mesh_3.h>
#include <CGAL/disable_warnings.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/optimize_mesh_3.h>

namespace parameters{

// -----------------------------------
// Perturb
// -----------------------------------

template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<internal::Perturb_options, internal_np::perturb_options_param_t> perturb(const CGAL_NP_CLASS& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    double time_limit = choose_parameter(get_parameter(np,internal_np::maximum_running_time),internal::undef_parameter);
    double sliver_bound = choose_parameter(get_parameter(np,internal_np::lower_sliver_bound),default_values_for_mesh_3::perturb_sliver_bound);

    internal::Perturb_options options(true);

    if ( internal::undef_parameter != time_limit)
        options.set_time_limit(time_limit);

    options.set_bound(sliver_bound);
    typedef Named_function_parameters<internal::Perturb_options, internal_np::perturb_options_param_t> Param;
    return Param(options);
}

template<typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
Named_function_parameters<internal::Perturb_options, internal_np::perturb_options_param_t> perturb(const CGAL_NP_CLASS& ... nps)
{
    return perturb(internal_np::combine_named_parameters(nps...));
}


inline Named_function_parameters<internal::Perturb_options, internal_np::perturb_options_param_t> no_perturb() {

    typedef Named_function_parameters<internal::Perturb_options, internal_np::perturb_options_param_t> Param;
    return Param(internal::Perturb_options(false));
}

// -----------------------------------
// Exude
// -----------------------------------
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<internal::Exude_options, internal_np::exude_options_param_t> exude(const CGAL_NP_CLASS& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    double time_limit = choose_parameter(get_parameter(np,internal_np::maximum_running_time),internal::undef_parameter);
    double sliver_bound = choose_parameter(get_parameter(np,internal_np::lower_sliver_bound),default_values_for_mesh_3::perturb_sliver_bound);

    internal::Exude_options options(true);

    if ( internal::undef_parameter != time_limit)
        options.set_time_limit(time_limit);
    options.set_bound(sliver_bound);
    typedef Named_function_parameters<internal::Exude_options, internal_np::exude_options_param_t> Param;

    return Param(options);
}

template<typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
Named_function_parameters<internal::Exude_options, internal_np::exude_options_param_t> exude(const CGAL_NP_CLASS& ... nps)
{
    return exude(internal_np::combine_named_parameters(nps...));
}

inline Named_function_parameters<internal::Exude_options, internal_np::exude_options_param_t> no_exude() {
    typedef Named_function_parameters<internal::Exude_options, internal_np::exude_options_param_t> Param;
    return Param(internal::Exude_options(false));
}

// -----------------------------------
// Odt
// -----------------------------------
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<internal::Odt_options, internal_np::odt_options_param_t> odt(const CGAL_NP_CLASS& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    double time_limit = choose_parameter(get_parameter(np,internal_np::maximum_running_time),0);
    double freeze_bound = choose_parameter(get_parameter(np,internal_np::vertex_freeze_bound),default_values_for_mesh_3::odt_freeze_ratio);
    double convergence = choose_parameter(get_parameter(np,internal_np::convergence_ratio), default_values_for_mesh_3::odt_convergence_ratio);
    int max_iteration_number = choose_parameter(get_parameter(np,internal_np::number_of_iterations), 0);
    internal::Odt_options options(true);

    options.set_time_limit(time_limit);
    options.set_bound(freeze_bound);
    options.set_convergence(convergence);
    options.set_max_iteration_number(max_iteration_number);
    typedef Named_function_parameters<internal::Odt_options, internal_np::odt_options_param_t> Param;
    return Param(options);
}

template<typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
Named_function_parameters<internal::Odt_options, internal_np::odt_options_param_t> odt(const CGAL_NP_CLASS& ... nps)
{
    return odt(internal_np::combine_named_parameters(nps...));
}

inline Named_function_parameters<internal::Odt_options, internal_np::odt_options_param_t> no_odt() {
    typedef Named_function_parameters<internal::Odt_options, internal_np::odt_options_param_t> Param;
    return Param(internal::Odt_options(false));
}

// -----------------------------------
// Lloyd
// -----------------------------------
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<internal::Lloyd_options, internal_np::lloyd_options_param_t> lloyd(const CGAL_NP_CLASS& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    double time_limit = choose_parameter(get_parameter(np,internal_np::maximum_running_time),0);
    double freeze_bound = choose_parameter(get_parameter(np,internal_np::vertex_freeze_bound),default_values_for_mesh_3::lloyd_freeze_ratio);
    double convergence = choose_parameter(get_parameter(np,internal_np::convergence_ratio), default_values_for_mesh_3::lloyd_convergence_ratio);
    int max_iteration_number = choose_parameter(get_parameter(np,internal_np::number_of_iterations), 0);
    internal::Lloyd_options options(true);

    options.set_time_limit(time_limit);
    options.set_bound(freeze_bound);
    options.set_convergence(convergence);
    options.set_max_iteration_number(max_iteration_number);

    typedef Named_function_parameters<internal::Lloyd_options, internal_np::lloyd_options_param_t> Param;
    return Param(options);
}


template<typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
Named_function_parameters<internal::Lloyd_options, internal_np::lloyd_options_param_t> lloyd(const CGAL_NP_CLASS& ... nps)
{
    return lloyd(internal_np::combine_named_parameters(nps...));
}

inline Named_function_parameters<internal::Lloyd_options, internal_np::lloyd_options_param_t> no_lloyd() {
    typedef Named_function_parameters<internal::Lloyd_options, internal_np::lloyd_options_param_t> Param;
    return Param(internal::Lloyd_options(false));
}

// -----------------------------------
// Manifold options ------------------
// -----------------------------------
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> manifold_options(const CGAL_NP_CLASS& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    int mesh_topology = choose_parameter(get_parameter(np, internal_np::mesh_topology_number), -1);
    internal::Manifold_options options;
    options.mesh_topology = mesh_topology;

    typedef Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> Param;
    return Param(options);
}


template<typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> manifold_options(const CGAL_NP_CLASS& ... nps)
{
    return manifold_options(internal_np::combine_named_parameters(nps...));
}

inline Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> manifold()
{
    typedef Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> Param;
  return Param(internal::Manifold_options(
          internal::Manifold_options::MANIFOLD));
}
inline Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> manifold_with_boundary()
{
    typedef Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> Param;
  return Param(internal::Manifold_options(
          internal::Manifold_options::MANIFOLD_WITH_BOUNDARY));
}
inline Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> non_manifold()
{
    typedef Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> Param;
  return Param(internal::Manifold_options(
          internal::Manifold_options::NON_MANIFOLD));
}

// -----------------------------------
// Mesh options
// -----------------------------------

// Undocumented Boost parameter for refine_mesh_3 and make_mesh_3.
// Allows to dump the mesh at given stage of the mesh generation
// algorithm.

template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<internal::Mesh_3_options, internal_np::mesh_param_t> mesh_3_options(const CGAL_NP_CLASS& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    internal::Mesh_3_options options;

    options.dump_after_init_prefix=choose_parameter(get_parameter(np, internal_np::dump_after_init_prefix_param), "");
    options.dump_after_refine_surface_prefix=choose_parameter(get_parameter(np, internal_np::dump_after_refine_surface_prefix_param), "");
    options.dump_after_refine_prefix=choose_parameter(get_parameter(np, internal_np::dump_after_refine_prefix_param), "");
    options.dump_after_glob_opt_prefix=choose_parameter(get_parameter(np, internal_np::dump_after_glob_opt_prefix_param), "");
    options.dump_after_perturb_prefix=choose_parameter(get_parameter(np, internal_np::dump_after_perturb_prefix_param), "");
    options.dump_after_exude_prefix=choose_parameter(get_parameter(np, internal_np::dump_after_refine_surface_prefix_param), "");
    options.number_of_initial_points=choose_parameter(get_parameter(np, internal_np::number_of_initial_points_param), -1);
    options.nonlinear_growth_of_balls = choose_parameter(get_parameter(np, internal_np::nonlinear_growth_of_balls_param), false);
    options.maximal_number_of_vertices=choose_parameter(get_parameter(np, internal_np::maximal_number_of_vertices_param), 0);
    options.pointer_to_error_code=choose_parameter(get_parameter(np, internal_np::pointer_to_error_code_param), ((Mesh_error_code*)0));
#ifndef CGAL_NO_ATOMIC
    options.pointer_to_stop_atomic_boolean=choose_parameter(get_parameter(np, internal_np::pointer_to_stop_atomic_boolean_param),
                                                            ((internal::Mesh_3_options::Pointer_to_stop_atomic_boolean_t)0));
#endif

    typedef Named_function_parameters<internal::Mesh_3_options, internal_np::mesh_param_t> Param;
    return Param(options);
}

template<typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
Named_function_parameters<internal::Mesh_3_options, internal_np::mesh_param_t> mesh_3_options(const CGAL_NP_CLASS& ... nps)
{
    return mesh_3_options(internal_np::combine_named_parameters(nps...));
}


} //namespace parameters

#endif //CGAL_PARAMETERS_H