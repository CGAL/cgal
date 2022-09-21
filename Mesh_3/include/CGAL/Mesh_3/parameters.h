// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_MESH_3_PARAMETERS_H
#define CGAL_MESH_3_PARAMETERS_H

#include <CGAL/license/Mesh_3.h>
#include <CGAL/Mesh_error_code.h>
#include <CGAL/Mesh_3/parameters_defaults.h>
#include <CGAL/Named_function_parameters.h>

namespace CGAL {

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

// Features
struct Features_options
{
  Features_options(bool b) : b_(b) {}
  bool features() const { return b_; }
private:
  bool b_;
};

// -----------------------------------
// Features generator
// -----------------------------------
// struct Features_option_generator
template <typename HasFeatures>
struct Features_options_generator {};

template<>
struct Features_options_generator<CGAL::Tag_true>
{
  Features_options operator()() { return Features_options(true); }
};

template<>
struct Features_options_generator<CGAL::Tag_false>
{
  Features_options operator()() { return Features_options(false); }
};

// struct Domain_features_generator is designed to handle cases where
// MeshDomain::Has_features is not a valid type
template< typename MeshDomain, bool MeshDomainHasHasFeatures >
struct Domain_features_generator {};

template< typename MeshDomain >
struct Domain_features_generator< MeshDomain, false >
{
  Features_options operator()()
  {
    return Features_options_generator<CGAL::Tag_false>()();
  }
};

template< typename MeshDomain >
struct Domain_features_generator< MeshDomain, true >
{
  Features_options operator()()
  {
    return Features_options_generator<typename MeshDomain::Has_features>()();
  }
};

} // end namespace internal



// -----------------------------------
// Reset_c3t3 (undocumented)
// -----------------------------------
inline
Named_function_parameters<bool, internal_np::do_reset_c3t3_t>
reset_c3t3()
{
  return Named_function_parameters<bool, internal_np::do_reset_c3t3_t>(true);
}

inline
Named_function_parameters<bool, internal_np::do_reset_c3t3_t>
no_reset_c3t3()
{
  return Named_function_parameters<bool, internal_np::do_reset_c3t3_t>(false);
}

// -----------------------------------
// Perturb
// -----------------------------------
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<internal::Perturb_options, internal_np::perturb_options_param_t>
perturb(const CGAL_NP_CLASS& np = parameters::default_values())
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
Named_function_parameters<internal::Perturb_options, internal_np::perturb_options_param_t>
perturb(const CGAL_NP_CLASS& ... nps)
{
  return perturb(internal_np::combine_named_parameters(nps...));
}


inline Named_function_parameters<internal::Perturb_options, internal_np::perturb_options_param_t>
no_perturb()
{
  typedef Named_function_parameters<internal::Perturb_options, internal_np::perturb_options_param_t> Param;
  return Param(internal::Perturb_options(false));
}

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_DEPRECATED
inline
Named_function_parameters<internal::Perturb_options, internal_np::perturb_options_param_t>
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
Named_function_parameters<internal::Exude_options, internal_np::exude_options_param_t>
exude(const CGAL_NP_CLASS& np = parameters::default_values())
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
Named_function_parameters<internal::Exude_options, internal_np::exude_options_param_t>
exude(const CGAL_NP_CLASS& ... nps)
{
  return exude(internal_np::combine_named_parameters(nps...));
}

inline Named_function_parameters<internal::Exude_options, internal_np::exude_options_param_t>
no_exude()
{
  typedef Named_function_parameters<internal::Exude_options, internal_np::exude_options_param_t> Param;
  return Param(internal::Exude_options(false));
}

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_DEPRECATED
inline
Named_function_parameters<internal::Exude_options, internal_np::exude_options_param_t>
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
Named_function_parameters<internal::Odt_options, internal_np::odt_options_param_t>
odt(const CGAL_NP_CLASS& np = parameters::default_values())
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
Named_function_parameters<internal::Odt_options, internal_np::odt_options_param_t>
odt(const CGAL_NP_CLASS& ... nps)
{
  return odt(internal_np::combine_named_parameters(nps...));
}

inline Named_function_parameters<internal::Odt_options, internal_np::odt_options_param_t>
no_odt()
{
  typedef Named_function_parameters<internal::Odt_options, internal_np::odt_options_param_t> Param;
  return Param(internal::Odt_options(false));
}

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_DEPRECATED
inline
Named_function_parameters<internal::Odt_options, internal_np::odt_options_param_t>
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
Named_function_parameters<internal::Lloyd_options, internal_np::lloyd_options_param_t>
lloyd(const CGAL_NP_CLASS& np = parameters::default_values())
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
Named_function_parameters<internal::Lloyd_options, internal_np::lloyd_options_param_t>
lloyd(const CGAL_NP_CLASS& ... nps)
{
  return lloyd(internal_np::combine_named_parameters(nps...));
}

inline Named_function_parameters<internal::Lloyd_options, internal_np::lloyd_options_param_t>
no_lloyd()
{
  typedef Named_function_parameters<internal::Lloyd_options, internal_np::lloyd_options_param_t> Param;
  return Param(internal::Lloyd_options(false));
}

#ifndef CGAL_NO_DEPRECATED_CODE
CGAL_DEPRECATED
inline
Named_function_parameters<internal::Lloyd_options, internal_np::lloyd_options_param_t>
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
Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t>
manifold_options(const CGAL_NP_CLASS& np = parameters::default_values())
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
Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t>
manifold_options(const CGAL_NP_CLASS& ... nps)
{
  return manifold_options(internal_np::combine_named_parameters(nps...));
}

inline Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t>
manifold()
{
  typedef Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> Param;
  return Param(internal::Manifold_options(internal::Manifold_options::MANIFOLD));
}
inline Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t>
manifold_with_boundary()
{
    typedef Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> Param;
  return Param(internal::Manifold_options(
          internal::Manifold_options::MANIFOLD_WITH_BOUNDARY));
}
inline Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t>
non_manifold()
{
  typedef Named_function_parameters<internal::Manifold_options, internal_np::manifold_param_t> Param;
  return Param(internal::Manifold_options(internal::Manifold_options::NON_MANIFOLD));
}

// -----------------------------------
// Mesh options
// -----------------------------------

// Undocumented parameter for refine_mesh_3 and make_mesh_3.
// Allows to dump the mesh at given stage of the mesh generation
// algorithm.
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Named_function_parameters<internal::Mesh_3_options, internal_np::mesh_param_t>
mesh_3_options(const CGAL_NP_CLASS& np = parameters::default_values())
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
Named_function_parameters<internal::Mesh_3_options, internal_np::mesh_param_t>
mesh_3_options(const CGAL_NP_CLASS& ... nps)
{
  return mesh_3_options(internal_np::combine_named_parameters(nps...));
}

// Undocumented parameter for refine_mesh_3 and make_mesh_3.
// Default Mesh_3_options: dump at every stage of the mesh generation.
inline
Named_function_parameters<internal::Mesh_3_options, internal_np::mesh_param_t>
mesh_3_dump()
{
  typedef Named_function_parameters<internal::Mesh_3_options, internal_np::mesh_param_t> Param;
  internal::Mesh_3_options options;

  options.dump_after_init_prefix = "mesh_dump_after_init";
  options.dump_after_refine_surface_prefix = "mesh_dump_after_refine_surface";
  options.dump_after_refine_prefix = "mesh_dump_after_refine";
  options.dump_after_glob_opt_prefix = "mesh_dump_after_glob_opt";
  options.dump_after_perturb_prefix = "mesh_dump_after_perturb";
  options.dump_after_exude_prefix = "mesh_dump_after_exude";

  return Param(options);
}

// -----------------------------------
// Features_options
// -----------------------------------
inline Named_function_parameters<internal::Features_options, internal_np::features_option_param_t>
features() {
  typedef Named_function_parameters<internal::Features_options, internal_np::features_option_param_t> Param;
  return Param(internal::Features_options(true));
}

inline Named_function_parameters<internal::Features_options, internal_np::features_option_param_t>
no_features() {
  typedef Named_function_parameters<internal::Features_options, internal_np::features_option_param_t> Param;
  return Param(internal::Features_options(false)); }

template < typename MeshDomain >
inline Named_function_parameters<internal::Features_options, internal_np::features_option_param_t>
features(const MeshDomain& /*domain*/)
{
typedef typename internal::Domain_features_generator<
  MeshDomain,
  CGAL::Mesh_3::internal::has_Has_features<MeshDomain>::value > Generator;

typedef Named_function_parameters<internal::Features_options, internal_np::features_option_param_t> Param;
return Param(Generator()());
}

} } //namespace CGAL::parameters

#endif //CGAL_MESH_3_PARAMETERS_H