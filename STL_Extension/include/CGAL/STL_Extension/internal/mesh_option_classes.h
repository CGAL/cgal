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

#include <functional>

#include <CGAL/STL_Extension/internal/Has_features.h>
#include <boost/iterator/function_output_iterator.hpp>

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
//   This process has two parameters : `initial_points_generator` and `initial_points`.
//   These two parameters are packed into a `Initial_points_generator_options`
//   that do not know the parameters types.
//   To remove the type of the `initial_points_generator` functor, two `std::function` are used.
//   To remove the type of the `initial_points` container, two `Input_const_iterator_interface` are used.

// A common interface for an iterator to a const value.
template <typename Value>
class Input_const_iterator_interface
{
public:
  virtual ~Input_const_iterator_interface() {}
  virtual const Value& operator*() = 0;
  virtual Input_const_iterator_interface<Value>* operator++() = 0;
  virtual bool operator!=(const Input_const_iterator_interface<Value>* other) const = 0;
  virtual Input_const_iterator_interface<Value>* clone() = 0;
};

// An iterator container that implements the `Input_const_iterator_interface` interface.
template <typename Value, typename Iterator>
struct Input_const_iterator_container
  : Input_const_iterator_interface<Value>
{
  typedef Input_const_iterator_container<Value, Iterator> Self;
public:
  Input_const_iterator_container(const Iterator& it) : it_(it) {}

  ~Input_const_iterator_container() override {}

  const Value& operator*() override
  {
    return *it_;
  }

  Input_const_iterator_interface<Value>* operator++() override
  {
    ++it_;
    return this;
  }

  bool operator!=(const Input_const_iterator_interface<Value>* other) const override
  {
    const Self* other_casted = dynamic_cast<const Self*>(other);
    if (other_casted == nullptr)
      return true;
    return it_ != other_casted->it_;
  }

  Input_const_iterator_interface<Value>* clone() override
  {
    return new Input_const_iterator_container<Value, Iterator>(it_);
  }

private:
  Iterator it_;
};

// Holds the two parameters `initial_points_generator` and `initial_points`,
// without knowing their types, into a single generator.
template <typename MeshDomain, typename C3t3>
struct Initial_points_generator_options
{
  typedef typename C3t3::Triangulation::Geom_traits::Weighted_point_3 Weighted_point_3;
  typedef typename MeshDomain::Index Index;
  typedef typename std::tuple<Weighted_point_3, int, Index> Value_type;
  typedef typename std::back_insert_iterator<std::vector<Value_type>> OutputIterator;

  template <typename Initial_points_generator, typename Initial_points>
  Initial_points_generator_options(const Initial_points_generator& generator, const Initial_points& initial_points, bool is_default = false)
    : initial_points_generator_no_number_of_points_(generator)
    , initial_points_generator_(generator)
    , is_default_(is_default && initial_points.size() == 0)
  {
    if (initial_points.size() == 0)
    {
      begin_it = nullptr;
      end_it = nullptr;
    }
    else
    {
      using Iterator_type = typename Initial_points::const_iterator;
      begin_it = new Input_const_iterator_container<Value_type, Iterator_type>(initial_points.cbegin());
      end_it = new Input_const_iterator_container<Value_type, Iterator_type>(initial_points.cend());
    }
  }

  ~Initial_points_generator_options()
  {
    if (begin_it != nullptr)
      delete begin_it;
    if (end_it != nullptr)
      delete end_it;
  }

  // Firstly, the `initial_points` are inserted, then, the `initial_points_generator` is called.
  OutputIterator operator()(OutputIterator pts, const MeshDomain& domain, const C3t3& c3t3) const
  {
    add_initial_points(pts);
    return initial_points_generator_no_number_of_points_(pts, domain, c3t3);
  }

  OutputIterator operator()(OutputIterator pts, const MeshDomain& domain, const C3t3& c3t3, int n) const
  {
    add_initial_points(pts);
    return initial_points_generator_(pts, domain, c3t3, n);
  }

  OutputIterator add_initial_points(OutputIterator pts) const
  {
    if (begin_it != nullptr && end_it != nullptr)
    {
      Input_const_iterator_interface<Value_type>* it_ptr = begin_it->clone();
      Input_const_iterator_interface<Value_type>& it = *(it_ptr);
      for (; it != end_it; ++it)
        *pts++ = *it;
      delete it_ptr;
    }
    return pts;
  }

  bool is_default() const { return is_default_; }

private:
  // The two functions holds the `initial_points_generator` functor
  const std::function<OutputIterator(OutputIterator&,const MeshDomain&,const C3t3&)> initial_points_generator_no_number_of_points_;
  const std::function<OutputIterator(OutputIterator&,const MeshDomain&,const C3t3&,int)> initial_points_generator_;
  // The two iterators holds the `initial_points` container
  Input_const_iterator_interface<Value_type>* begin_it;
  Input_const_iterator_interface<Value_type>* end_it;
  // Is the option a default-constructed one ?
  const bool is_default_;
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

// Evaluate the two parameters `initial_points_generator` and `initial_points`
// and returns the appropriate `Initial_points_generator_options`.
// If no generator and no initial points, then use the domain's construct_initial_points_object.
template <typename MeshDomain, typename C3t3>
struct Initial_points_generator_generator
{
  typedef typename CGAL::parameters::internal::Initial_points_generator_options<MeshDomain, C3t3> Initial_points_generator_options;
  typedef typename Initial_points_generator_options::Value_type     Value_type;
  typedef typename Initial_points_generator_options::OutputIterator OutputIterator;

  struct Initial_points_generator_domain_traductor
  {
    OutputIterator operator()(OutputIterator pts, const MeshDomain& domain, const C3t3& c3t3)
    {
      // Use boost to easily create an output iterator.
      // This iterator take the domain's construct_initial_points_object output : an std::pair<Point_3, Index>
      // and outputs an std::tuple<Weighted_point_3, dimension, Index>
      // As points are on the surfaces by construction, dimension is always 2.
      typename C3t3::Triangulation::Geom_traits::Construct_weighted_point_3 cwp =
          c3t3.triangulation().geom_traits().construct_weighted_point_3_object();
      domain.construct_initial_points_object()(
             boost::make_function_output_iterator([&](const auto& domain_generated_point) {
               *pts++ = std::make_tuple(cwp(domain_generated_point.first), 2, domain_generated_point.second);
             }));
      return pts;
    }
    OutputIterator operator()(OutputIterator pts, const MeshDomain& domain, const C3t3& c3t3, int n)
    {
      typename C3t3::Triangulation::Geom_traits::Construct_weighted_point_3 cwp =
          c3t3.triangulation().geom_traits().construct_weighted_point_3_object();
      domain.construct_initial_points_object()(
             boost::make_function_output_iterator([&](const auto& domain_generated_point) {
               *pts++ = std::make_tuple(cwp(domain_generated_point.first), 2, domain_generated_point.second);
             }), n);
      return pts;
    }
  };

  struct Initial_points_generator_empty
  {
    OutputIterator operator()(OutputIterator pts, const MeshDomain& , const C3t3& )
    { return pts; }
    OutputIterator operator()(OutputIterator pts, const MeshDomain& , const C3t3& , int )
    { return pts; }
  };

  // With a custom InitialPointsGenerator
  template <typename InitialPointsGenerator, typename InitalPointsRange>
  Initial_points_generator_options operator()(const InitialPointsGenerator& initial_points_generator, const InitalPointsRange& input_features)
  {
    return Initial_points_generator_options(initial_points_generator, input_features, false);
  }

  template <typename InitialPointsGenerator>
  Initial_points_generator_options operator()(const InitialPointsGenerator& initial_points_generator)
  {
    std::vector<Value_type> empty_input_features;
    return operator()(initial_points_generator, empty_input_features);
  }

  // Without a custom InitialPointsGenerator
  template <typename InitalPointsRange>
  Initial_points_generator_options operator()(const Null_functor&, const InitalPointsRange& input_features)
  {
    // The domain's construct_initial_points_object is called only if input_features is empty
    if (input_features.size() == 0) {
      return Initial_points_generator_options(Initial_points_generator_domain_traductor(), input_features, true);
    }
    return Initial_points_generator_options(Initial_points_generator_empty(), input_features, true);
  }

  // Default construction
  Initial_points_generator_options operator()()
  {
    return operator()(Null_functor());
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
