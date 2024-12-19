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
#include <CGAL/type_traits.h>
#include <CGAL/STL_Extension/internal/tuple_like_helpers.h>

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
// Holds the two parameters `initial_points_generator` and `initial_points`,
// without knowing their types, into a single generator.
template <typename MeshDomain,
          typename C3t3,
          typename InitialPointsRange = CGAL::Default
>
struct Initialization_options
{
  using Point = typename C3t3::Triangulation::Point;
  using Default_initial_point_type
          = std::tuple<Point, int, typename MeshDomain::Index>;
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

  struct Output_function_ref {
    // This reference-like type uses type erasure to store a reference to a callable
    //
    // See the video "Breaking Dependencies - C++ Type Erasure - The Implementation Details"
    //   by Klaus Iglberger at CppCon 2022, from time code 49:57.
    // URL: https://youtu.be/qn6OqefuH08?si=YzhwpgNLur8_jOeC&t=2997"

    using Erased_call_function_pointer_type = void(*)(void*, const Default_initial_point_type&);

    // store the address of the callable
    void* const f_ = nullptr;
    // and the call function (the non-capturing lambda generated by the templated constructor)
    Erased_call_function_pointer_type const call_ = nullptr;

    template <typename Function,
              typename = std::enable_if_t<!std::is_same_v<CGAL::cpp20::remove_cvref_t<Function>,
                                                          Output_function_ref>
                                          >
              >
    Output_function_ref(Function&& f)
      : f_(std::addressof(f))
      , call_( [](void* f, const Default_initial_point_type& p) {
                 using F = CGAL::cpp20::remove_cvref_t<Function>;
                 auto* real_f = static_cast<F*>(f);
                 (*real_f)(p);
               } )
    {
    }

    template <typename Tuple_like>
    void operator()(Tuple_like&& p) const
    {
      using Tuple_like_no_cvref = CGAL::cpp20::remove_cvref_t<Tuple_like>;
      if constexpr (CGAL::STL_Extension::internal::tuple_like_of_size_2<Tuple_like_no_cvref>) {
        const auto& [pt, index] = p;
        call_(f_, Default_initial_point_type(pt, 2, index));
      } else if constexpr (std::is_same_v<Tuple_like_no_cvref, Default_initial_point_type>) {
        call_(f_, std::forward<Tuple_like>(p));
      } else {
        const auto& [pt, dim, index] = p;
        call_(f_, Default_initial_point_type(pt, dim, index));
      }
    }
  }; // end of struct Output_function_ref

  using Point_output_function_iterator = boost::function_output_iterator<Output_function_ref>;

  struct Generator_ref { // type-erased reference to a generator, same as Output_function_ref
    using Erased_call_function_pointer_type = Point_output_function_iterator(*)(void*, Point_output_function_iterator, const int);

    void * const generator_ = nullptr;
    Erased_call_function_pointer_type const call_ = nullptr;

    template <typename Generator,
              typename = std::enable_if_t<!std::is_same_v<CGAL::cpp20::remove_cvref_t<Generator>,
                                                          Generator_ref>
                                          >
              >
    Generator_ref(Generator&& generator)
      : generator_(std::addressof(generator))
      , call_( [](void* g, Point_output_function_iterator oit, const int n) {
                 using Real_generator = CGAL::cpp20::remove_cvref_t<Generator>;
                 auto* real_g = static_cast<Real_generator*>(g);
                 return (*real_g)(oit, n);
               } )
    {
    }

    Generator_ref() = default;

    Point_output_function_iterator operator()(Point_output_function_iterator oit, const int n) const
    {
      return call_(generator_, oit, n);
    }

    Point_output_function_iterator operator()(Point_output_function_iterator oit, const int n)
    {
      return call_(generator_, oit, n);
    }

    bool operator==(std::nullptr_t) const { return generator_ == nullptr; }
  }; // end of struct Generator_ref

  Initialization_options()
  {}

  template <typename Initial_points_generator>
  Initialization_options(Initial_points_generator& generator,
                         const Initial_points_range& initial_points)
    : initial_points_generator_(std::forward<Initial_points_generator>(generator))
    , begin_it(cbegin(initial_points))
    , end_it(cend(initial_points))
  {}

  template <typename Self, typename OutputIterator>
  static OutputIterator call_operator(Self& self, OutputIterator pts_it, const int n)
  {
    // add initial_points
    pts_it = std::copy(self.begin_it, self.end_it, pts_it);

    if(self.initial_points_generator_ == nullptr) {
      return pts_it;
    }

    // Now, create an output iterator type-erasing the type of `pts_it`...
    auto output_to_pts_it = [&](const Default_initial_point_type& p) { *pts_it++ = p; };
    Output_function_ref function_ref{ output_to_pts_it }; // maintains a non-const reference to pts_it
    Point_output_function_iterator output_iterator{ function_ref };

    // ...and use the type-erased output iterator with the type-erased generator.
    self.initial_points_generator_(output_iterator, n);
    return pts_it;
  }

  template<typename OutputIterator>
  OutputIterator operator()(OutputIterator pts, const int n = 0)
  {
    return call_operator(*this, pts, n);
  }

  template<typename OutputIterator>
  OutputIterator operator()(OutputIterator pts, const int n = 0) const
  {
    return call_operator(*this, pts, n);
  }

  bool is_default() const
  {
    return begin_it == end_it && initial_points_generator_ == nullptr;
  }

private:
  Generator_ref initial_points_generator_; //reference that type-erases the generator type

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
