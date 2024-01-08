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

#ifndef CGAL_NO_STATIC_ASSERTION_TESTS
#include <CGAL/basic.h>
#endif

#include <CGAL/tags.h>
#include <CGAL/STL_Extension/internal/mesh_option_classes.h>

#include <boost/mpl/has_xxx.hpp>

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

enum all_default_t { all_default };

// define enum types and values for new named parameters
#define CGAL_add_named_parameter(X, Y, Z)            \
  enum X { Y };
#define CGAL_add_named_parameter_with_compatibility(X, Y, Z)            \
  enum X { Y };
#define CGAL_add_named_parameter_with_compatibility_cref_only(X, Y, Z)            \
  enum X { Y };
#define CGAL_add_named_parameter_with_compatibility_ref_only(X, Y, Z)            \
  enum X { Y };
#define CGAL_add_extra_named_parameter_with_compatibility(X, Y, Z)
#include <CGAL/STL_Extension/internal/parameters_interface.h>
#undef CGAL_add_named_parameter
#undef CGAL_add_named_parameter_with_compatibility
#undef CGAL_add_named_parameter_with_compatibility_cref_only
#undef CGAL_add_named_parameter_with_compatibility_ref_only
#undef CGAL_add_extra_named_parameter_with_compatibility

template <typename T, typename Tag, typename Base>
struct Named_params_impl : Base
{
  typename std::conditional<std::is_copy_constructible<T>::value,
                            T, std::reference_wrapper<const T> >::type v; // copy of the parameter if copyable
  Named_params_impl(const T& v, const Base& b)
    : Base(b)
    , v(v)
  {}
};

// partial specialization for base class of the recursive nesting
template <typename T, typename Tag>
struct Named_params_impl<T, Tag, No_property>
{
  typename std::conditional<std::is_copy_constructible<T>::value,
                            T, std::reference_wrapper<const T> >::type v; // copy of the parameter if copyable
  Named_params_impl(const T& v)
    : v(v)
  {}
};

// Helper class to get the type of a named parameter pack given a query tag
template <typename NP, typename Query_tag>
struct Get_param;

template< typename T, typename Tag, typename Query_tag>
struct Get_param< Named_params_impl<T, Tag, No_property>, Query_tag >
{
  typedef Param_not_found type;
  typedef Param_not_found reference;
};

template< typename T, typename Tag, typename Base>
struct Get_param< Named_params_impl<T, Tag, Base>, Tag >
{
  typedef typename std::conditional<std::is_copy_constructible<T>::value,
                                    T, std::reference_wrapper<const T> >::type type;
  typedef typename std::conditional<std::is_copy_constructible<T>::value,
                                    T, const T&>::type reference;
};

template< typename T, typename Tag>
struct Get_param< Named_params_impl<T, Tag, No_property>, Tag >
{
  typedef typename std::conditional<std::is_copy_constructible<T>::value,
                                    T, std::reference_wrapper<const T> >::type type;
  typedef typename std::conditional<std::is_copy_constructible<T>::value,
                                    T, const T&>::type reference;
};

template< typename T, typename Tag, typename Base>
struct Get_param< Named_params_impl<std::reference_wrapper<T>, Tag, Base>, Tag >
{
  typedef std::reference_wrapper<T> type;
  typedef T& reference;
};

template< typename T, typename Tag>
struct Get_param< Named_params_impl<std::reference_wrapper<T>, Tag, No_property>, Tag >
{
  typedef std::reference_wrapper<T> type;
  typedef T& reference;
};


template< typename T, typename Tag, typename Base, typename Query_tag>
struct Get_param< Named_params_impl<T,Tag,Base>, Query_tag>
{
  typedef typename Get_param<typename Base::base, Query_tag>::type type;
  typedef typename Get_param<typename Base::base, Query_tag>::reference reference;
};

// helper to choose the default
template <typename Query_tag, typename NP, typename D>
struct Lookup_named_param_def
{
  typedef typename internal_np::Get_param<typename NP::base, Query_tag>::type NP_type;
  typedef typename internal_np::Get_param<typename NP::base, Query_tag>::reference NP_reference;

  typedef std::conditional_t<
    std::is_same_v<NP_type, internal_np::Param_not_found>,
    D, NP_type>
  type;

  typedef std::conditional_t<
    std::is_same_v<NP_reference, internal_np::Param_not_found>,
    D&, NP_reference>
  reference;
};

// helper function to extract the value from a named parameter pack given a query tag
template <typename T, typename Tag, typename Base>
typename std::conditional<std::is_copy_constructible<T>::value,
                          T, std::reference_wrapper<const T> >::type
get_parameter_impl(const Named_params_impl<T, Tag, Base>& np, Tag)
{
  return np.v;
}

template< typename T, typename Tag, typename Query_tag>
Param_not_found get_parameter_impl(const Named_params_impl<T, Tag, No_property>&, Query_tag)
{
  return Param_not_found();
}

template< typename T, typename Tag>
typename std::conditional<std::is_copy_constructible<T>::value,
                          T, std::reference_wrapper<const T> >::type
get_parameter_impl(const Named_params_impl<T, Tag, No_property>& np, Tag)
{
  return np.v;
}

template <typename T, typename Tag, typename Base, typename Query_tag>
typename Get_param<Named_params_impl<T, Tag, Base>, Query_tag>::type
get_parameter_impl(const Named_params_impl<T, Tag, Base>& np, Query_tag tag)
{
#ifndef CGAL_NO_STATIC_ASSERTION_TEST
  static_assert(!std::is_same<Query_tag, Tag>::value);
#endif
  return get_parameter_impl(static_cast<const typename Base::base&>(np), tag);
}


// helper for getting references
template <class T>
const T& get_reference(const T& t)
{
  return t;
}

template <class T>
T& get_reference(const std::reference_wrapper<T>& r)
{
  return r.get();
}

// helper function to extract the reference from a named parameter pack given a query tag
template <typename T, typename Tag, typename Base>
typename std::conditional<std::is_copy_constructible<T>::value,
                          T, const T& >::type
get_parameter_reference_impl(const Named_params_impl<T, Tag, Base>& np, Tag)
{
  return get_reference(np.v);
}

template< typename T, typename Tag, typename Query_tag>
Param_not_found
get_parameter_reference_impl(const Named_params_impl<T, Tag, No_property>&, Query_tag)
{
  return Param_not_found();
}

template< typename T, typename Tag>
typename std::conditional<std::is_copy_constructible<T>::value,
                          T, const T& >::type
get_parameter_reference_impl(const Named_params_impl<T, Tag, No_property>& np, Tag)
{
  return get_reference(np.v);
}

template <typename T, typename Tag, typename Base>
T&
get_parameter_reference_impl(const Named_params_impl<std::reference_wrapper<T>, Tag, Base>& np, Tag)
{
  return np.v.get();
}

template< typename T, typename Tag>
T&
get_parameter_reference_impl(const Named_params_impl<std::reference_wrapper<T>, Tag, No_property>& np, Tag)
{
  return np.v.get();
}

template <typename T, typename Tag, typename Base, typename Query_tag>
typename Get_param<Named_params_impl<T, Tag, Base>, Query_tag>::reference
get_parameter_reference_impl(const Named_params_impl<T, Tag, Base>& np, Query_tag tag)
{
  static_assert(!std::is_same<Query_tag, Tag>::value);
  return get_parameter_reference_impl(static_cast<const typename Base::base&>(np), tag);
}

} // end of internal_np namespace

template <typename T, typename Tag, typename Base = internal_np::No_property>
struct Named_function_parameters;

namespace parameters{

typedef Named_function_parameters<bool, internal_np::all_default_t>  Default_named_parameters;

Default_named_parameters
inline default_values();

// function to extract a parameter
template <typename T, typename Tag, typename Base, typename Query_tag>
typename internal_np::Get_param<internal_np::Named_params_impl<T, Tag, Base>, Query_tag>::type
get_parameter(const Named_function_parameters<T, Tag, Base>& np, Query_tag tag)
{
  return internal_np::get_parameter_impl(static_cast<const internal_np::Named_params_impl<T, Tag, Base>&>(np), tag);
}

template <typename T, typename Tag, typename Base, typename Query_tag>
typename internal_np::Get_param<internal_np::Named_params_impl<T, Tag, Base>, Query_tag>::reference
get_parameter_reference(const Named_function_parameters<T, Tag, Base>& np, Query_tag tag)
{
  return internal_np::get_parameter_reference_impl(
    static_cast<const internal_np::Named_params_impl<T, Tag, Base>&>(np),
    tag);
}

// Two parameters, non-trivial default value
template <typename D>
D& choose_parameter(const internal_np::Param_not_found&, D& d)
{
  return d;
}

template <typename D>
const D& choose_parameter(const internal_np::Param_not_found&, const D& d)
{
  return d;
}

template <typename D>
D choose_parameter(const internal_np::Param_not_found&, D&& d)
{
  return std::forward<D>(d);
}

template <typename T, typename D>
T& choose_parameter(T& t, D&)
{
  return t;
}

template <typename T, typename D>
const T& choose_parameter(const T& t, const D&)
{
  return t;
}

// single parameter so that we can avoid a default construction
template <typename D>
D choose_parameter(const internal_np::Param_not_found&)
{
  return D();
}

template <typename D, typename T>
const T& choose_parameter(const T& t)
{
  return t;
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

  Named_function_parameters() : base(T()) {}
  Named_function_parameters(const T& v) : base(v) {}
  Named_function_parameters(const T& v, const Base& b) : base(v, b) {}

// create the functions for new named parameters and the one imported boost
// used to concatenate several parameters
#define CGAL_add_named_parameter(X, Y, Z)                             \
  template<typename K>                                                \
  Named_function_parameters<K, internal_np::X, self>                  \
  Z(const K& k) const                                                 \
  {                                                                   \
    typedef Named_function_parameters<K, internal_np::X, self> Params;\
    return Params(k, *this);                                          \
  }
#define CGAL_add_named_parameter_with_compatibility(X, Y, Z)          \
  template<typename K>                                                \
  Named_function_parameters<K, internal_np::X, self>                  \
  Z(const K& k) const                                                 \
  {                                                                   \
    typedef Named_function_parameters<K, internal_np::X, self> Params;\
    return Params(k, *this);                                          \
  }
#define CGAL_add_named_parameter_with_compatibility_cref_only(X, Y, Z) \
  template<typename K>                                                \
  Named_function_parameters<std::reference_wrapper<const K>,          \
                            internal_np::X, self>                     \
  Z(const K& k) const                                                 \
  {                                                                   \
    typedef Named_function_parameters<std::reference_wrapper<const K>,\
                                      internal_np::X, self> Params;   \
    return Params(std::cref(k), *this);                               \
  }
#define CGAL_add_named_parameter_with_compatibility_ref_only(X, Y, Z) \
  template<typename K>                                                \
  Named_function_parameters<std::reference_wrapper<K>,                \
                            internal_np::X, self>                     \
  Z(K& k) const                                                       \
  {                                                                   \
    typedef Named_function_parameters<std::reference_wrapper<K>,      \
                                      internal_np::X, self> Params;   \
    return Params(std::ref(k), *this);                                \
  }
#define CGAL_add_extra_named_parameter_with_compatibility(X, Y, Z)    \
  template<typename K>                                                \
  Named_function_parameters<K, internal_np::X, self>                  \
  Z(const K& k) const                                                 \
  {                                                                   \
    typedef Named_function_parameters<K, internal_np::X, self> Params;\
    return Params(k, *this);                                          \
  }
#include <CGAL/STL_Extension/internal/parameters_interface.h>
#undef CGAL_add_named_parameter
#undef CGAL_add_named_parameter_with_compatibility
#undef CGAL_add_named_parameter_with_compatibility_cref_only
#undef CGAL_add_named_parameter_with_compatibility_ref_only
#undef CGAL_add_extra_named_parameter_with_compatibility

// inject mesh specific named parameter functions
#define CGAL_NP_BASE self
#define CGAL_NP_BUILD(P, V) P(V, *this)

#include <CGAL/STL_Extension/internal/mesh_parameters_interface.h>

#undef CGAL_NP_BASE
#undef CGAL_NP_BUILD

  template <typename OT, typename OTag>
  Named_function_parameters<OT, OTag, self>
  combine(const Named_function_parameters<OT,OTag>& np) const
  {
    return Named_function_parameters<OT, OTag, self>(np.v,*this);
  }

  template <typename OT, typename OTag, typename ... NPS>
  auto
  combine(const Named_function_parameters<OT,OTag>& np, const NPS& ... nps) const
  {
    return Named_function_parameters<OT, OTag, self>(np.v,*this).combine(nps...);
  }

  // typedef for SFINAE
  typedef int CGAL_Named_function_parameters_class;
};

namespace parameters {

Default_named_parameters
inline default_values()
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
  Named_function_parameters<K, Tag>
  operator()(const K& p) const
  {
    typedef Named_function_parameters<K, Tag> Params;
    return Params(p);
  }

  template <typename K>
  Named_function_parameters<K, Tag>
  operator=(const K& p) const
  {
    typedef Named_function_parameters<K, Tag> Params;
    return Params(p);
  }
};

template <class Tag>
struct Boost_parameter_compatibility_wrapper<Tag, true, true>
{
  template <typename K>
  Named_function_parameters<std::reference_wrapper<const K>, Tag>
  operator()(const K& p) const
  {
    typedef Named_function_parameters<std::reference_wrapper<const K>, Tag> Params;
    return Params(std::cref(p));
  }

  template <typename K>
  Named_function_parameters<std::reference_wrapper<const K>, Tag>
  operator=(const K& p) const
  {
    typedef Named_function_parameters<std::reference_wrapper<const K>, Tag> Params;
    return Params(std::cref(p));
  }
};

template <class Tag>
struct Boost_parameter_compatibility_wrapper<Tag, true, false>
{
  template <typename K>
  Named_function_parameters<std::reference_wrapper<K>, Tag>
  operator()(K& p) const
  {
    typedef Named_function_parameters<std::reference_wrapper<K>, Tag> Params;
    return Params(std::ref(p));
  }

  template <typename K>
  Named_function_parameters<std::reference_wrapper<K>, Tag>
  operator=(std::reference_wrapper<K> p) const
  {
    typedef Named_function_parameters<std::reference_wrapper<K>, Tag> Params;
    return Params(std::ref(p));
  }
};

// define free functions and Boost_parameter_compatibility_wrapper for named parameters
#define CGAL_add_named_parameter(X, Y, Z)        \
  template <typename K>                        \
  Named_function_parameters<K, internal_np::X>                  \
  Z(const K& p)                                \
  {                                            \
    typedef Named_function_parameters<K, internal_np::X> Params;\
    return Params(p);                          \
  }

#define CGAL_add_named_parameter_with_compatibility(X, Y, Z)        \
  const Boost_parameter_compatibility_wrapper<internal_np::X> Z;
#define CGAL_add_named_parameter_with_compatibility_cref_only(X, Y, Z)        \
  const Boost_parameter_compatibility_wrapper<internal_np::X, true, true> Z;
#define CGAL_add_named_parameter_with_compatibility_ref_only(X, Y, Z)        \
  const Boost_parameter_compatibility_wrapper<internal_np::X, true, false> Z;
#define CGAL_add_extra_named_parameter_with_compatibility(X, Y, Z)        \
  const Boost_parameter_compatibility_wrapper<internal_np::X> Z;
#include <CGAL/STL_Extension/internal/parameters_interface.h>
#undef CGAL_add_named_parameter
#undef CGAL_add_extra_named_parameter_with_compatibility
#undef CGAL_add_named_parameter_with_compatibility
#undef CGAL_add_named_parameter_with_compatibility_cref_only
#undef CGAL_add_named_parameter_with_compatibility_ref_only

// Version with three parameters for dynamic property maps
template <typename D, typename Dynamic_tag, typename PolygonMesh>
D choose_parameter(const internal_np::Param_not_found&, Dynamic_tag tag, PolygonMesh& pm)
{
  return get(tag, pm);
}

template <typename D, typename T, typename Dynamic_tag, typename PolygonMesh>
const T& choose_parameter(const T& t, Dynamic_tag, PolygonMesh&)
{
  return t;
}

template <class NamedParameters, class Parameter>
struct is_default_parameter
{
  typedef typename internal_np::Lookup_named_param_def<Parameter,
                                                       NamedParameters,
                                                       internal_np::Param_not_found>::type NP_type;

  static const bool value = std::is_same<NP_type, internal_np::Param_not_found>::value;

  typedef CGAL::Boolean_tag<value> type;
};


// code used to make sure all options passed are used by a function
namespace authorized_parameters_impl
{

template <class ... Tag>
struct Tag_wrapper{};

template <class TagAllowed, class ... TagsAllowed, class Tag>
constexpr
bool is_tag_present(Tag_wrapper<TagAllowed, TagsAllowed...>, Tag)
{
  if (std::is_same_v<TagAllowed, Tag>)
    return true;
  else
    return is_tag_present(Tag_wrapper<TagsAllowed...>(), Tag());
}

template <class Tag>
constexpr
bool is_tag_present(Tag_wrapper<Tag>, Tag)
{
  return true;
}

template <class TagAllowed, class Tag>
constexpr
bool is_tag_present(Tag_wrapper<TagAllowed>, Tag)
{
  return false;
}

template <class ... TagsAllowed, class T, class Tag>
constexpr
bool authorized_options_rec(const Named_function_parameters<T, Tag>&)
{
  return is_tag_present(Tag_wrapper<TagsAllowed...>(), Tag());
}

template <class ... TagsAllowed, class T, class Tag, class Base>
constexpr
bool authorized_options_rec(const Named_function_parameters<T, Tag, Base>& np)
{
  if (is_tag_present(Tag_wrapper<TagsAllowed...>(), Tag()))
    return authorized_options_rec<TagsAllowed...>(static_cast<const Base&>(np));
  return false;
}

}// impl namespace

template <class ... TagsAllowed, class Named_function_parameters>
constexpr
bool authorized_options(const Named_function_parameters& np)
{
#ifndef CGAL_DISABLE_NAMED_FUNCTION_PARAMETERS_CHECKS
  return authorized_parameters_impl::authorized_options_rec
    <internal_np::all_default_t, TagsAllowed...>(np);
#else
  return true;
#endif
}

} // end of parameters namespace

#ifndef CGAL_NO_DEPRECATED_CODE
namespace Polygon_mesh_processing {

namespace parameters = CGAL::parameters;

}
#endif

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

// For disambiguation using SFINAE
BOOST_MPL_HAS_XXX_TRAIT_DEF(CGAL_Named_function_parameters_class)
template<class T>
inline constexpr bool is_named_function_parameter = has_CGAL_Named_function_parameters_class<T>::value;

#endif // CGAL_BOOST_FUNCTION_PARAMS_HPP
