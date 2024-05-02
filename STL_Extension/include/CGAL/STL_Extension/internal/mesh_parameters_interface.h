// Copyright (c) 2017  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//

// List of named parameters special functions used in meshing packages of CGAL
// no guard on purpose as they are injected both in the parameter namespace
// and the Named_function_parameters class.


// -----------------------------------
// Reset_c3t3 (undocumented)
// -----------------------------------
inline
Named_function_parameters<bool, ::CGAL::internal_np::do_reset_c3t3_t, CGAL_NP_BASE>
reset_c3t3()
{
  typedef Named_function_parameters<bool, ::CGAL::internal_np::do_reset_c3t3_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param, true);
}

inline
Named_function_parameters<bool, ::CGAL::internal_np::do_reset_c3t3_t, CGAL_NP_BASE>
no_reset_c3t3()
{
  typedef Named_function_parameters<bool, ::CGAL::internal_np::do_reset_c3t3_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param, false);
}

// -----------------------------------
// Perturb
// -----------------------------------
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<::CGAL::parameters::internal::Perturb_options, ::CGAL::internal_np::perturb_options_param_t, CGAL_NP_BASE>
perturb(const CGAL_NP_CLASS& np = parameters::default_values())
{
  using ::CGAL::parameters::choose_parameter;
  using ::CGAL::parameters::get_parameter;
  double time_limit = choose_parameter(get_parameter(np,::CGAL::internal_np::maximum_running_time),::CGAL::parameters::internal::undef_parameter);
  double sliver_bound = choose_parameter(get_parameter(np,::CGAL::internal_np::lower_sliver_bound),::CGAL::parameters::default_values_for_mesh_3::perturb_sliver_bound);

  ::CGAL::parameters::internal::Perturb_options options(true);

  if ( ::CGAL::parameters::internal::undef_parameter != time_limit)
      options.set_time_limit(time_limit);

  options.set_bound(sliver_bound);
  typedef Named_function_parameters<::CGAL::parameters::internal::Perturb_options, ::CGAL::internal_np::perturb_options_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param, options);
}

template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1, typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2, typename ... NP>
Named_function_parameters<::CGAL::parameters::internal::Perturb_options, ::CGAL::internal_np::perturb_options_param_t, CGAL_NP_BASE>
perturb(const CGAL_NP_CLASS_1&  np1, const CGAL_NP_CLASS_2&  np2, const NP& ... nps)
{
  return perturb(::CGAL::internal_np::combine_named_parameters(np1, np2, nps...));
}


inline Named_function_parameters<::CGAL::parameters::internal::Perturb_options, ::CGAL::internal_np::perturb_options_param_t, CGAL_NP_BASE>
no_perturb()
{
  typedef Named_function_parameters<::CGAL::parameters::internal::Perturb_options, ::CGAL::internal_np::perturb_options_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param,::CGAL::parameters::internal::Perturb_options(false));
}

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_DEPRECATED
inline
Named_function_parameters<::CGAL::parameters::internal::Perturb_options, ::CGAL::internal_np::perturb_options_param_t, CGAL_NP_BASE>
perturb(double time_limit_,
        double sliver_bound_=0)
{
  return perturb(time_limit(time_limit_).
                 sliver_bound(sliver_bound_));
}
#endif

// -----------------------------------
// Exude
// -----------------------------------
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<::CGAL::parameters::internal::Exude_options, ::CGAL::internal_np::exude_options_param_t, CGAL_NP_BASE>
exude(const CGAL_NP_CLASS& np = parameters::default_values())
{
  using ::CGAL::parameters::choose_parameter;
  using ::CGAL::parameters::get_parameter;
  double time_limit = choose_parameter(get_parameter(np,::CGAL::internal_np::maximum_running_time),::CGAL::parameters::internal::undef_parameter);
  double sliver_bound = choose_parameter(get_parameter(np,::CGAL::internal_np::lower_sliver_bound),::CGAL::parameters::default_values_for_mesh_3::exude_sliver_bound);

  ::CGAL::parameters::internal::Exude_options options(true);

  if ( ::CGAL::parameters::internal::undef_parameter != time_limit)
      options.set_time_limit(time_limit);
  options.set_bound(sliver_bound);
  typedef Named_function_parameters<::CGAL::parameters::internal::Exude_options, ::CGAL::internal_np::exude_options_param_t, CGAL_NP_BASE> Param;

  return CGAL_NP_BUILD(Param, options);
}

template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1, typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2, typename ... NP>
Named_function_parameters<::CGAL::parameters::internal::Exude_options, ::CGAL::internal_np::exude_options_param_t, CGAL_NP_BASE>
exude(const CGAL_NP_CLASS_1&  np1, const CGAL_NP_CLASS_2&  np2, const NP& ... nps)
{
  return exude(::CGAL::internal_np::combine_named_parameters(np1, np2, nps...));
}

inline Named_function_parameters<::CGAL::parameters::internal::Exude_options, ::CGAL::internal_np::exude_options_param_t, CGAL_NP_BASE>
no_exude()
{
  typedef Named_function_parameters<::CGAL::parameters::internal::Exude_options, ::CGAL::internal_np::exude_options_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param,::CGAL::parameters::internal::Exude_options(false));
}

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_DEPRECATED
inline
Named_function_parameters<::CGAL::parameters::internal::Exude_options, ::CGAL::internal_np::exude_options_param_t, CGAL_NP_BASE>
exude(double time_limit_,
      double sliver_bound_ = 0)
{
  return exude(time_limit(time_limit_).sliver_bound(sliver_bound_));
}
#endif

// -----------------------------------
// Odt
// -----------------------------------
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<::CGAL::parameters::internal::Odt_options, ::CGAL::internal_np::odt_options_param_t, CGAL_NP_BASE>
odt(const CGAL_NP_CLASS& np = parameters::default_values())
{
  using ::CGAL::parameters::choose_parameter;
  using ::CGAL::parameters::get_parameter;
  double time_limit = choose_parameter(get_parameter(np,::CGAL::internal_np::maximum_running_time),0);
  double freeze_bound = choose_parameter(get_parameter(np,::CGAL::internal_np::vertex_freeze_bound),::CGAL::parameters::default_values_for_mesh_3::odt_freeze_ratio);
  double convergence = choose_parameter(get_parameter(np,::CGAL::internal_np::convergence_ratio), ::CGAL::parameters::default_values_for_mesh_3::odt_convergence_ratio);
  std::size_t max_iteration_number = choose_parameter(get_parameter(np,::CGAL::internal_np::number_of_iterations), 0);
  ::CGAL::parameters::internal::Odt_options options(true);

  options.set_time_limit(time_limit);
  options.set_bound(freeze_bound);
  options.set_convergence(convergence);
  options.set_max_iteration_number(max_iteration_number);
  typedef Named_function_parameters<::CGAL::parameters::internal::Odt_options, ::CGAL::internal_np::odt_options_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param,options);
}

template<typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
Named_function_parameters<::CGAL::parameters::internal::Odt_options, ::CGAL::internal_np::odt_options_param_t, CGAL_NP_BASE>
odt(const CGAL_NP_CLASS& ... nps)
{
  return odt(::CGAL::internal_np::combine_named_parameters(nps...));
}

inline Named_function_parameters<::CGAL::parameters::internal::Odt_options, ::CGAL::internal_np::odt_options_param_t, CGAL_NP_BASE>
no_odt()
{
  typedef Named_function_parameters<::CGAL::parameters::internal::Odt_options, ::CGAL::internal_np::odt_options_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param,::CGAL::parameters::internal::Odt_options(false));
}

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_DEPRECATED
inline
Named_function_parameters<::CGAL::parameters::internal::Odt_options, ::CGAL::internal_np::odt_options_param_t, CGAL_NP_BASE>
odt(double time_limit_,
    std::size_t max_iteration_number_ = 0,
    double convergence_ = 0.02,
    double freeze_bound_ = 0.01,
    bool do_freeze_ = true)
{
  return odt(time_limit(time_limit_).
             max_iteration_number(max_iteration_number_).
             convergence(convergence_).
             freeze_bound(freeze_bound_).
             do_freeze(do_freeze_));
}
#endif

// -----------------------------------
// Lloyd
// -----------------------------------
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<::CGAL::parameters::internal::Lloyd_options, ::CGAL::internal_np::lloyd_options_param_t, CGAL_NP_BASE>
lloyd(const CGAL_NP_CLASS& np = parameters::default_values())
{
  using ::CGAL::parameters::choose_parameter;
  using ::CGAL::parameters::get_parameter;
  double time_limit = choose_parameter(get_parameter(np,::CGAL::internal_np::maximum_running_time),0);
  double freeze_bound = choose_parameter(get_parameter(np,::CGAL::internal_np::vertex_freeze_bound),::CGAL::parameters::default_values_for_mesh_3::lloyd_freeze_ratio);
  double convergence = choose_parameter(get_parameter(np,::CGAL::internal_np::convergence_ratio), ::CGAL::parameters::default_values_for_mesh_3::lloyd_convergence_ratio);
  std::size_t max_iteration_number = choose_parameter(get_parameter(np,::CGAL::internal_np::number_of_iterations), 0);
  ::CGAL::parameters::internal::Lloyd_options options(true);

  options.set_time_limit(time_limit);
  options.set_bound(freeze_bound);
  options.set_convergence(convergence);
  options.set_max_iteration_number(max_iteration_number);

  typedef Named_function_parameters<::CGAL::parameters::internal::Lloyd_options, ::CGAL::internal_np::lloyd_options_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param, options);
}

template<typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
Named_function_parameters<::CGAL::parameters::internal::Lloyd_options, ::CGAL::internal_np::lloyd_options_param_t, CGAL_NP_BASE>
lloyd(const CGAL_NP_CLASS& ... nps)
{
  return lloyd(::CGAL::internal_np::combine_named_parameters(nps...));
}

inline Named_function_parameters<::CGAL::parameters::internal::Lloyd_options, ::CGAL::internal_np::lloyd_options_param_t, CGAL_NP_BASE>
no_lloyd()
{
  typedef Named_function_parameters<::CGAL::parameters::internal::Lloyd_options, ::CGAL::internal_np::lloyd_options_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param, ::CGAL::parameters::internal::Lloyd_options(false));
}

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_DEPRECATED
inline
Named_function_parameters<::CGAL::parameters::internal::Lloyd_options, ::CGAL::internal_np::lloyd_options_param_t, CGAL_NP_BASE>
lloyd(double time_limit_,
      std::size_t max_iteration_number_ = 0,
      double convergence_ = 0.02,
      double freeze_bound_ = 0.01,
      bool do_freeze_= true)
{
  return lloyd(time_limit(time_limit_).
               max_iteration_number(max_iteration_number_).
               convergence(convergence_).
               freeze_bound(freeze_bound_).
               do_freeze(do_freeze_));
}
#endif

// -----------------------------------
// Manifold options
// -----------------------------------
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<::CGAL::parameters::internal::Manifold_options, ::CGAL::internal_np::manifold_param_t, CGAL_NP_BASE>
manifold_options(const CGAL_NP_CLASS& np = parameters::default_values())
{
  using ::CGAL::parameters::choose_parameter;
  using ::CGAL::parameters::get_parameter;
  int mesh_topology = choose_parameter(get_parameter(np, ::CGAL::internal_np::mesh_topology_number), -1);
  ::CGAL::parameters::internal::Manifold_options options;
  options.mesh_topology = mesh_topology;

  typedef Named_function_parameters<::CGAL::parameters::internal::Manifold_options, ::CGAL::internal_np::manifold_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param, options);
}


template<typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
Named_function_parameters<::CGAL::parameters::internal::Manifold_options, ::CGAL::internal_np::manifold_param_t, CGAL_NP_BASE>
manifold_options(const CGAL_NP_CLASS& ... nps)
{
  return manifold_options(::CGAL::internal_np::combine_named_parameters(nps...));
}

inline Named_function_parameters<::CGAL::parameters::internal::Manifold_options, ::CGAL::internal_np::manifold_param_t, CGAL_NP_BASE>
manifold()
{
  typedef Named_function_parameters<::CGAL::parameters::internal::Manifold_options, ::CGAL::internal_np::manifold_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param, ::CGAL::parameters::internal::Manifold_options(::CGAL::parameters::internal::Manifold_options::MANIFOLD));
}
inline Named_function_parameters<::CGAL::parameters::internal::Manifold_options, ::CGAL::internal_np::manifold_param_t, CGAL_NP_BASE>
manifold_with_boundary()
{
  typedef Named_function_parameters<::CGAL::parameters::internal::Manifold_options, ::CGAL::internal_np::manifold_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param,::CGAL::parameters::internal::Manifold_options(
                               ::CGAL::parameters::internal::Manifold_options::MANIFOLD_WITH_BOUNDARY));
}
inline Named_function_parameters<::CGAL::parameters::internal::Manifold_options, ::CGAL::internal_np::manifold_param_t, CGAL_NP_BASE>
non_manifold()
{
  typedef Named_function_parameters<::CGAL::parameters::internal::Manifold_options, ::CGAL::internal_np::manifold_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param, ::CGAL::parameters::internal::Manifold_options(::CGAL::parameters::internal::Manifold_options::NON_MANIFOLD));
}

// -----------------------------------
// Mesh options
// -----------------------------------

// Undocumented parameter for refine_mesh_3 and make_mesh_3.
// Allows to dump the mesh at given stage of the mesh generation
// algorithm.
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<::CGAL::parameters::internal::Mesh_3_options, ::CGAL::internal_np::mesh_param_t, CGAL_NP_BASE>
mesh_3_options(const CGAL_NP_CLASS& np = parameters::default_values())
{
  using ::CGAL::parameters::choose_parameter;
  using ::CGAL::parameters::get_parameter;
  ::CGAL::parameters::internal::Mesh_3_options options;

  options.dump_after_init_prefix=choose_parameter(get_parameter(np, ::CGAL::internal_np::dump_after_init_prefix_param), "");
  options.dump_after_refine_surface_prefix=choose_parameter(get_parameter(np, ::CGAL::internal_np::dump_after_refine_surface_prefix_param), "");
  options.dump_after_refine_prefix=choose_parameter(get_parameter(np, ::CGAL::internal_np::dump_after_refine_prefix_param), "");
  options.dump_after_glob_opt_prefix=choose_parameter(get_parameter(np, ::CGAL::internal_np::dump_after_glob_opt_prefix_param), "");
  options.dump_after_perturb_prefix=choose_parameter(get_parameter(np, ::CGAL::internal_np::dump_after_perturb_prefix_param), "");
  options.dump_after_exude_prefix=choose_parameter(get_parameter(np, ::CGAL::internal_np::dump_after_exude_prefix_param), "");
  options.number_of_initial_points=choose_parameter(get_parameter(np, ::CGAL::internal_np::number_of_initial_points_param), -1);
  options.nonlinear_growth_of_balls = choose_parameter(get_parameter(np, ::CGAL::internal_np::nonlinear_growth_of_balls_param), false);
  options.maximal_number_of_vertices=choose_parameter(get_parameter(np, ::CGAL::internal_np::maximal_number_of_vertices_param), 0);
  options.pointer_to_error_code=choose_parameter(get_parameter(np, ::CGAL::internal_np::pointer_to_error_code_param), ((Mesh_error_code*)0));
#ifndef CGAL_NO_ATOMIC
  options.pointer_to_stop_atomic_boolean=choose_parameter(get_parameter(np, ::CGAL::internal_np::pointer_to_stop_atomic_boolean_param),
                                                          ((::CGAL::parameters::internal::Mesh_3_options::Pointer_to_stop_atomic_boolean_t)0));
#endif

  typedef Named_function_parameters<::CGAL::parameters::internal::Mesh_3_options, ::CGAL::internal_np::mesh_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param, options);
}

template<typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
Named_function_parameters<::CGAL::parameters::internal::Mesh_3_options, ::CGAL::internal_np::mesh_param_t, CGAL_NP_BASE>
mesh_3_options(const CGAL_NP_CLASS& ... nps)
{
  return mesh_3_options(::CGAL::internal_np::combine_named_parameters(nps...));
}

// Undocumented parameter for refine_mesh_3 and make_mesh_3.
// Default Mesh_3_options: dump at every stage of the mesh generation.
inline
Named_function_parameters<::CGAL::parameters::internal::Mesh_3_options, ::CGAL::internal_np::mesh_param_t, CGAL_NP_BASE>
mesh_3_dump()
{
  typedef Named_function_parameters<::CGAL::parameters::internal::Mesh_3_options, ::CGAL::internal_np::mesh_param_t, CGAL_NP_BASE> Param;
  ::CGAL::parameters::internal::Mesh_3_options options;

  options.dump_after_init_prefix = "mesh_dump_after_init";
  options.dump_after_refine_surface_prefix = "mesh_dump_after_refine_surface";
  options.dump_after_refine_prefix = "mesh_dump_after_refine";
  options.dump_after_glob_opt_prefix = "mesh_dump_after_glob_opt";
  options.dump_after_perturb_prefix = "mesh_dump_after_perturb";
  options.dump_after_exude_prefix = "mesh_dump_after_exude";

  return CGAL_NP_BUILD(Param, options);
}

// -----------------------------------
// Features_options
// -----------------------------------
inline Named_function_parameters<::CGAL::parameters::internal::Features_options, ::CGAL::internal_np::features_option_param_t, CGAL_NP_BASE>
features() {
  typedef Named_function_parameters<::CGAL::parameters::internal::Features_options, ::CGAL::internal_np::features_option_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param, ::CGAL::parameters::internal::Features_options(true));
}

inline Named_function_parameters<::CGAL::parameters::internal::Features_options, ::CGAL::internal_np::features_option_param_t, CGAL_NP_BASE>
no_features() {
  typedef Named_function_parameters<::CGAL::parameters::internal::Features_options, ::CGAL::internal_np::features_option_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param, ::CGAL::parameters::internal::Features_options(false)); }

template < typename MeshDomain >
inline Named_function_parameters<::CGAL::parameters::internal::Features_options, ::CGAL::internal_np::features_option_param_t, CGAL_NP_BASE>
features(const MeshDomain& /*domain*/)
{
  typedef typename ::CGAL::parameters::internal::Domain_features_generator<MeshDomain,
    CGAL::internal::has_Has_features<MeshDomain>::value > Generator;

  typedef Named_function_parameters<::CGAL::parameters::internal::Features_options, ::CGAL::internal_np::features_option_param_t, CGAL_NP_BASE> Param;
  return CGAL_NP_BUILD(Param,Generator()());
}
