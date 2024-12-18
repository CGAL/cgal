// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//

#ifndef CGAL_MESH_OPTION_CLASSES_H
#define CGAL_MESH_OPTION_CLASSES_H

#include <CGAL/STL_Extension/internal/Has_features.h>
#include <CGAL/Default.h>

#include <boost/iterator/function_output_iterator.hpp>

#include <iterator>
#include <optional>
#include <type_traits>

namespace CGAL {

enum Mesh_error_code {
  CGAL_MESH_3_NO_ERROR = 0,
  CGAL_MESH_3_MAXIMAL_NUMBER_OF_VERTICES_REACHED,
  CGAL_MESH_3_STOPPED
};

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

  bool is_max_iteration_number_set() const { return max_it_nb_ != std::size_t(undef_parameter); }
  void set_max_iteration_number(std::size_t i) { max_it_nb_ = i; }
  std::size_t max_iteration_number() const { return max_it_nb_; }

private:
  double convergence_;
  std::size_t max_it_nb_;
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

// Mesh initialization
struct Dummy_initial_points_generator
{
  template<typename OutputIterator>
  OutputIterator operator()(OutputIterator oit, const int) const { return oit; }
};

// Holds the two parameters `initial_points_generator` and `initial_points`,
// without knowing their types, into a single generator.
template <typename MeshDomain,
          typename C3t3,
          typename InitialPointsGenerator = CGAL::Default,
          typename InitialPointsRange = CGAL::Default
>
struct Initialization_options
{
  using Default_generator = Dummy_initial_points_generator;
  using Initial_points_generator
          = typename CGAL::Default::Get<InitialPointsGenerator, Default_generator>::type;
  using Default_initial_point_type
          = std::tuple<typename C3t3::Triangulation::Point, int, typename MeshDomain::Index>;
  using Initial_points_range
          = typename CGAL::Default::Get<InitialPointsRange, std::vector<Default_initial_point_type>>::type;

  template <typename Range>
  static auto cbegin(Range&& range) {
    return std::cbegin(std::forward<Range>(range));
  }

  template <typename Range>
  static auto cend(Range&& range) {
    return std::cend(std::forward<Range>(range));
  }
  using Initial_points_const_iterator = decltype(cbegin(std::declval<Initial_points_range>()));

  Initialization_options()
  {}

  Initialization_options(const Initial_points_generator& generator,
                         const Initial_points_range& initial_points)
    : initial_points_generator_(generator)
    , begin_it(cbegin(initial_points))
    , end_it(cend(initial_points))
  {}

  template<typename OutputIterator>
  OutputIterator operator()(OutputIterator pts, const int n = 0) const
  {
    // add initial_points
    for (Initial_points_const_iterator it = begin_it; it != end_it; ++it)
      *pts++ = *it;

    return initial_points_generator_(pts, n);
  }

  bool is_default() const
  {
    return begin_it == end_it
        && std::is_same_v<InitialPointsGenerator, Default_generator>;
  }

private:
  Initial_points_generator initial_points_generator_;

  // The two iterators point to the `initial_points` container
  Initial_points_const_iterator begin_it;
  Initial_points_const_iterator end_it;
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


namespace default_values_for_mesh_3 {

const double time_limit = 0.;

// exude_mesh_3
const double exude_sliver_bound = 0.;

// perturb_mesh_3
const double perturb_sliver_bound = 0.;

// lloyd_optimize_mesh_3
const double lloyd_freeze_ratio = 0.01;
const double lloyd_convergence_ratio = 0.02;

// odt_optimize_mesh_3
const double odt_freeze_ratio = 0.01;
const double odt_convergence_ratio = 0.02;

// global optimizers
const bool do_freeze = true;

} } } //namespace CGAL::parameters::def

#endif
