// Copyright (c) 2019  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_NAMED_FUNCTION_PARAMETERS_H
#define CGAL_NAMED_FUNCTION_PARAMETERS_H

#include <CGAL/type_traits.h>
#ifndef CGAL_NO_STATIC_ASSERTION_TESTS
#include <CGAL/basic.h>
#endif

#include <CGAL/tags.h>
#include <CGAL/STL_Extension/internal/mesh_option_classes.h>

#include <boost/mpl/has_xxx.hpp>

#include <functional>
#include <type_traits>
#include <utility>

#define CGAL_NP_TEMPLATE_PARAMETERS NP_T=bool, typename NP_Tag=CGAL::internal_np::all_default_t, typename NP_Base=CGAL::internal_np::No_property
#define CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT NP_T, typename NP_Tag, typename NP_Base
#define CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1 NP_T1, typename NP_Tag1, typename NP_Base1
#define CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2 NP_T2, typename NP_Tag2, typename NP_Base2
#define CGAL_NP_CLASS CGAL::Named_function_parameters<NP_T,NP_Tag,NP_Base>

#define CGAL_NP_TEMPLATE_PARAMETERS_1 NP_T1=bool, typename NP_Tag1=CGAL::internal_np::all_default_t, typename NP_Base1=CGAL::internal_np::No_property
#define CGAL_NP_CLASS_1 CGAL::Named_function_parameters<NP_T1,NP_Tag1,NP_Base1>
#define CGAL_NP_TEMPLATE_PARAMETERS_2 NP_T2=bool, typename NP_Tag2=CGAL::internal_np::all_default_t, typename NP_Base2=CGAL::internal_np::No_property
#define CGAL_NP_CLASS_2 CGAL::Named_function_parameters<NP_T2,NP_Tag2,NP_Base2>

#define CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC NP_T, typename ... NP_Tag, typename ... NP_Base

namespace CGAL {
namespace internal_np{

struct No_property {};
struct Param_not_found {};

template <typename T>
inline constexpr bool is_param_not_found_v = std::is_same_v<CGAL::cpp20::remove_cvref_t<T>, Param_not_found>;

enum all_default_t { all_default };

// helper for getting references
template <class T>
T get_reference(const T& t)
{
  return t;
}

template <class T>
T& get_reference(const std::reference_wrapper<T>& r)
{
  return r.get();
}

// define enum types and values for new named parameters
#define CGAL_add_named_parameter(X, Y, Z) \
  enum X { Y };
#include <CGAL/STL_Extension/internal/parameters_interface.h>

} // end namespace internal_np

// forward-declaration of Named_function_parameters
template <typename T, typename Tag, typename Base = internal_np::No_property>
struct Named_function_parameters;

namespace internal_np {
template <typename T, typename Tag, typename Base>
struct Named_params_impl : Base
{
  typename std::conditional<std::is_copy_constructible<T>::value,
                            T, std::reference_wrapper<const T> >::type v; // copy of the parameter if copyable
  Named_params_impl(const T& v, const Base& b)
    : Base(b)
    , v(v)
  {}

  constexpr decltype(auto) parameter(Tag) const noexcept { return v; }
  static constexpr bool has_parameter(Tag) noexcept { return true; }
  using Base::parameter;
  using Base::has_parameter;
};

// partial specialization for base class of the recursive nesting
template <typename T, typename Tag>
struct Named_params_impl<T, Tag, No_property>
{
  typename std::conditional<std::is_copy_constructible<T>::value,
                            T, std::reference_wrapper<const T> >::type v; // copy of the parameter if copyable
  constexpr Named_params_impl(const T& v)
    : v(v)
  {}
  constexpr decltype(auto) parameter(Tag) const noexcept { return v; }
  static constexpr bool has_parameter(Tag) noexcept { return true; }
};

// Helper class to get the type of a named parameter pack given a query tag
template <typename NP, typename Query_tag>
struct Get_param;

template <typename T, typename Tag, typename Base, typename Query_tag>
struct Get_param<Named_params_impl<T, Tag, Base>, Query_tag>
{
  using type = decltype(std::declval<Named_function_parameters<T, Tag, Base>>().parameter(Query_tag{}));
  using reference =
      decltype(get_reference(std::declval<Named_function_parameters<T, Tag, Base>>().parameter(Query_tag{})));
};

// helper to choose the default
template <typename Query_tag, typename NP, typename D>
struct Lookup_named_param_def
{
  typedef typename internal_np::Get_param<typename NP::base, Query_tag>::type NP_type;
  typedef typename internal_np::Get_param<typename NP::base, Query_tag>::reference NP_reference;

  typedef std::conditional_t<
    internal_np::is_param_not_found_v<NP_type>,
    D, NP_type>
  type;

  typedef std::conditional_t<
    internal_np::is_param_not_found_v<NP_reference>,
    D&, NP_reference>
  reference;
};

} // end of internal_np namespace

template <typename U, typename T, typename Tag, typename Base, typename Query_tag>
constexpr decltype(auto) parameter_or([[maybe_unused]] const Named_function_parameters<T, Tag, Base>& np,
                                      [[maybe_unused]] Query_tag tag,
                                      [[maybe_unused]] U&& default_value)
{
  return np.parameter_or(tag, std::forward<U>(default_value));
}

template <typename U, typename T, typename Tag, typename Base, typename Query_tag>
constexpr decltype(auto) parameter_or([[maybe_unused]] const Named_function_parameters<T, Tag, Base>& np,
                                      [[maybe_unused]] Query_tag tag)
{
  return np.template parameter_or<U>(tag);
}

namespace parameters{

typedef Named_function_parameters<bool, internal_np::all_default_t>  Default_named_parameters;

inline constexpr Default_named_parameters default_values();

// function to extract a parameter
template <typename T, typename Tag, typename Base, typename Query_tag>
constexpr decltype(auto)
get_parameter(const Named_function_parameters<T, Tag, Base>& np, Query_tag tag)
{
  return np.parameter(tag);
}

template <typename T, typename Tag, typename Base, typename Query_tag>
constexpr decltype(auto)
get_parameter_reference(const Named_function_parameters<T, Tag, Base>& np, Query_tag tag)
{
  return internal_np::get_reference(np.parameter(tag));
}

// Two parameters, non-trivial default value
template <typename T, typename D>
constexpr decltype(auto) choose_parameter([[maybe_unused]] T&& t, [[maybe_unused]] D&& d) {
  if constexpr (internal_np::is_param_not_found_v<T>) {
    return std::forward<D>(d);
  } else {
    return std::forward<T>(t);
  }
}

// single parameter so that we can avoid a default construction
template <typename D, typename T>
constexpr decltype(auto) choose_parameter([[maybe_unused]] T&& t)
{
  if constexpr (internal_np::is_param_not_found_v<T>) {
    return D{};
  } else {
    return std::forward<T>(t);
  }
}

// version with a dynamic property tag with initialization
template <typename T, typename Tag, typename Graph, typename V>
constexpr decltype(auto) choose_parameter([[maybe_unused]] T&& t,
                                          [[maybe_unused]] Tag tag,
                                          [[maybe_unused]] Graph& graph,
                                          [[maybe_unused]] const V& default_value) {
  if constexpr (internal_np::is_param_not_found_v<T>) {
    return get(tag, graph, default_value);
  } else {
    return std::forward<T>(t);
  }
}

template <typename T, typename Tag, typename Graph>
constexpr decltype(auto)
choose_parameter([[maybe_unused]] T&& t, [[maybe_unused]] Tag tag, [[maybe_unused]] Graph& graph)
{
  if constexpr (internal_np::is_param_not_found_v<T>) {
    return get(tag, graph);
  } else {
    return std::forward<T>(t);
  }
}

} // parameters namespace

namespace internal_np {

template <typename Tag, typename K, typename ... NPS>
auto
combine_named_parameters(const Named_function_parameters<K, Tag>& np, const NPS& ... nps)
{
  return np.combine(nps ...);
}

} // end of internal_np namespace

template <typename T, typename Tag, typename Base>
struct Named_function_parameters
  : internal_np::Named_params_impl<T, Tag, Base>
{
  typedef internal_np::Named_params_impl<T, Tag, Base> base;
  typedef Named_function_parameters<T, Tag, Base> self;

  using base::parameter;
  using base::has_parameter;
  constexpr auto parameter(...) const { return internal_np::Param_not_found(); }
  static constexpr bool has_parameter(...) { return false; }

  template <typename D, typename Query_tag>
  constexpr decltype(auto) parameter_or([[maybe_unused]] Query_tag tag, [[maybe_unused]] D&& default_value) const {
    if constexpr (has_parameter(Query_tag())) {
      return parameter(tag);
    } else {
      return std::forward<D>(default_value);
    }
  }

  template <typename D, typename Query_tag>
  constexpr decltype(auto) parameter_or([[maybe_unused]] Query_tag tag) const {
    if constexpr (has_parameter(Query_tag())) {
      return parameter(tag);
    } else {
      return D{};
    }
  }

  constexpr Named_function_parameters() : base(T()) {}
  constexpr Named_function_parameters(const T& v) : base(v) {}
  constexpr Named_function_parameters(const T& v, const Base& b) : base(v, b) {}

// create the functions for new named parameters and the one imported boost
// used to concatenate several parameters
#define CGAL_add_named_parameter(X, Y, Z)                             \
  template<typename K>                                                \
  constexpr auto Z(const K& k) const                                  \
  {                                                                   \
    using Params = Named_function_parameters<K, internal_np::X, self>;\
    return Params(k, *this);                                          \
  }
#define CGAL_add_named_parameter_with_compatibility(X, Y, Z)          \
  CGAL_add_named_parameter(X, Y, Z)
#define CGAL_add_extra_named_parameter_with_compatibility(X, Y, Z)    \
  CGAL_add_named_parameter(X, Y, Z)
#define CGAL_add_named_parameter_with_compatibility_cref_only(X, Y, Z)\
  template<typename K>                                                \
  constexpr auto Z(const K& k) const                                  \
  {                                                                   \
    using Params =                                                    \
        Named_function_parameters<std::reference_wrapper<const K>,    \
                                  internal_np::X, self>;              \
    return Params(std::cref(k), *this);                               \
  }
#define CGAL_add_named_parameter_with_compatibility_ref_only(X, Y, Z) \
  template<typename K>                                                \
  constexpr auto Z(K& k) const                                        \
  {                                                                   \
    using Params =                                                    \
        Named_function_parameters<std::reference_wrapper<K>,          \
                                  internal_np::X, self>;              \
    return Params(std::ref(k), *this);                                \
  }
#include <CGAL/STL_Extension/internal/parameters_interface.h>

// inject mesh specific named parameter functions
#define CGAL_NP_BASE self
#define CGAL_NP_BUILD(P, V) P(V, *this)

#include <CGAL/STL_Extension/internal/mesh_parameters_interface.h>

#undef CGAL_NP_BASE
#undef CGAL_NP_BUILD

  template <typename OT, typename OTag>
  constexpr Named_function_parameters<OT, OTag, self>
  combine(const Named_function_parameters<OT,OTag>& np) const
  {
    return Named_function_parameters<OT, OTag, self>(np.v,*this);
  }

  template <typename OT, typename OTag, typename ... NPS>
  constexpr auto
  combine(const Named_function_parameters<OT,OTag>& np, const NPS& ... nps) const
  {
    return Named_function_parameters<OT, OTag, self>(np.v,*this).combine(nps...);
  }

  // typedef for SFINAE
  typedef int CGAL_Named_function_parameters_class;
};

namespace parameters {

inline constexpr Default_named_parameters default_values()
{
  return Default_named_parameters();
}

#ifndef CGAL_NO_DEPRECATED_CODE
Default_named_parameters
inline all_default()
{
  return Default_named_parameters();
}
#endif

template <class Tag, bool ref_only = false, bool ref_is_const = false>
struct Boost_parameter_compatibility_wrapper
{
  template <typename K>
  constexpr auto operator()(K&& p) const
  {
    using KK = cpp20::remove_cvref_t<cpp20::unwrap_ref_decay_t<K>>;
    if constexpr (ref_only)
    {
      if constexpr (ref_is_const)
      {
        using Params = Named_function_parameters<std::reference_wrapper<const KK>, Tag>;
        const auto& ref = cpp20::unwrap_reference_t<K>(p);
        return Params{std::cref(ref)};
      }
      else
      {
        using Params = Named_function_parameters<std::reference_wrapper<KK>, Tag>;
        auto& ref = cpp20::unwrap_reference_t<K>(p);
        return Params{std::ref(ref)};
      }
    }
    else
    {
      using Params = Named_function_parameters<K, Tag>;
      return Params{std::forward<K>(p)};
    }
  }

  template <typename K>
  constexpr auto operator=(K&& p) const
  {
    return operator()(std::forward<K>(p));
  }
};

// define free functions and Boost_parameter_compatibility_wrapper for named parameters
#define CGAL_add_named_parameter(X, Y, Z)                         \
  template <typename K>                                           \
  constexpr auto Z(const K& p) {                                  \
    using Params = Named_function_parameters<K, internal_np::X>;  \
    return Params{p};                                             \
  }

#define CGAL_add_named_parameter_with_compatibility(X, Y, Z)        \
  inline constexpr Boost_parameter_compatibility_wrapper<internal_np::X> Z;
#define CGAL_add_named_parameter_with_compatibility_cref_only(X, Y, Z)        \
  inline constexpr Boost_parameter_compatibility_wrapper<internal_np::X, true, true> Z;
#define CGAL_add_named_parameter_with_compatibility_ref_only(X, Y, Z)        \
  inline constexpr Boost_parameter_compatibility_wrapper<internal_np::X, true, false> Z;
#define CGAL_add_extra_named_parameter_with_compatibility(X, Y, Z)        \
  inline constexpr Boost_parameter_compatibility_wrapper<internal_np::X> Z;
#include <CGAL/STL_Extension/internal/parameters_interface.h>

template <class NamedParameters, class Parameter>
using is_default_parameter = Boolean_tag<!NamedParameters::has_parameter(Parameter{})>;

} // end of parameters namespace

#ifndef CGAL_NO_DEPRECATED_CODE
namespace Polygon_mesh_processing {

namespace parameters = CGAL::parameters;

}
#endif

// For disambiguation using SFINAE
BOOST_MPL_HAS_XXX_TRAIT_DEF(CGAL_Named_function_parameters_class)
template<class T>
inline constexpr bool is_named_function_parameter = has_CGAL_Named_function_parameters_class<T>::value;

} //namespace CGAL

#ifndef CGAL_NO_STATIC_ASSERTION_TESTS
// code added to avoid silent runtime issues in non-updated code
namespace boost
{
  template <typename T, typename Tag, typename Base, typename Tag2, bool B = false>
  void get_param(CGAL::Named_function_parameters<T,Tag,Base>, Tag2)
  {
    static_assert(B && "You must use CGAL::parameters::get_parameter instead of boost::get_param");
  }
}
#endif

#endif // CGAL_BOOST_FUNCTION_PARAMS_HPP
