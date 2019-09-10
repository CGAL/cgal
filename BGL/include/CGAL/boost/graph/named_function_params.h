//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// This file is part of the Boost Graph Library
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// $URL$
// $Id$
// SPDX-License-Identifier: BSL-1.0
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_BOOST_GRAPH_NAMED_FUNCTION_PARAMS_H
#define CGAL_BOOST_GRAPH_NAMED_FUNCTION_PARAMS_H

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>

#include <CGAL/boost/graph/properties.h>
#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4172) // returning address of local variable or temporary
#endif

#include <boost/graph/named_function_params.hpp>

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif


#include <boost/type_traits/is_same.hpp>
#include <boost/version.hpp>

// An explanation about the version hackery below: There is no real
// API to introduce custom properties to the Graph API and the
// internals have changed with Boost Version 1.51 and changes aren't
// backward compatible. To work around that we carry around two
// versions of cgal_bgl_named_params. One imitates the pre 1.51
// bgl_named_params, the newer one hooks into the API through
// inheritance and addition of the some partial specializations.

#if BOOST_VERSION < 105100
namespace boost{
  typedef detail::error_property_not_found param_not_found;

  template <typename Tag, typename Args, typename Def>
  struct lookup_named_param_def {
    typedef Def type;
    static const Def& get(const Args&, const Def& def) {return def;}
  };

  template <typename T, typename Tag, typename Base, typename Def>
  struct lookup_named_param_def<Tag, bgl_named_params<T, Tag, Base>, Def> {
    typedef T type;
    static const type& get(const bgl_named_params<T, Tag, Base>& p, const Def&) {
      return p.m_value;
    }
  };

  template <typename Tag1, typename T, typename Tag, typename Base, typename Def>
  struct lookup_named_param_def<Tag1, bgl_named_params<T, Tag, Base>, Def> {
    typedef typename lookup_named_param_def<Tag1, Base, Def>::type type;
    static const type& get(const bgl_named_params<T, Tag, Base>& p, const Def& def) {
      return lookup_named_param_def<Tag1, Base, Def>::get(p.m_base, def);
    }
  };
} //end of namespace boost
#endif

#define CGAL_BGL_NP_TEMPLATE_PARAMETERS T, typename Tag, typename Base
#define CGAL_BGL_NP_CLASS CGAL::cgal_bgl_named_params<T,Tag,Base>

namespace CGAL {
namespace internal_np{
enum all_default_t { all_default }; //cannot use macro because it takes no argument

// for uniformity we import them in this namespace. Note that
// it is an import so that if we use the named parameter function
// from boost it will work
using boost::vertex_index_t;
using boost::vertex_index;
using boost::graph_visitor_t;
using boost::visitor;

// define enum types and values for new named parameters
#define CGAL_add_named_parameter(X, Y, Z)            \
  enum X { Y };
#include <CGAL/boost/graph/parameters_interface.h>
#undef CGAL_add_named_parameter

}//internal_np

  template <typename T, typename Tag, typename Base = boost::no_property>
  struct cgal_bgl_named_params : boost::bgl_named_params<T, Tag, Base>
  {
    typedef boost::bgl_named_params<T, Tag, Base> base;
    typedef cgal_bgl_named_params self;

    cgal_bgl_named_params(T v = T()) : base(v) {}
    cgal_bgl_named_params(T v, const Base& b) : base(v, b) {}
    cgal_bgl_named_params<bool, internal_np::all_default_t, self>
    all_default() const
    {
      typedef cgal_bgl_named_params<bool, internal_np::all_default_t, self> Params;
      return Params(*this);
    }

// create the functions for new named parameters and the one imported boost
// used to concatenate several parameters
#define CGAL_add_named_parameter(X, Y, Z)                          \
  template<typename K>                                           \
  cgal_bgl_named_params<K, internal_np::X, self>                  \
  Z(const K& k) const                                            \
  {                                                              \
    typedef cgal_bgl_named_params<K, internal_np::X, self> Params;\
    return Params(k, *this);                                     \
  }
#include <CGAL/boost/graph/parameters_interface.h>
#include <CGAL/boost/graph/boost_parameters_interface.h>
#undef CGAL_add_named_parameter
  };

namespace parameters {

  cgal_bgl_named_params<bool, internal_np::all_default_t>
  inline all_default()
  {
    typedef cgal_bgl_named_params<bool, internal_np::all_default_t> Params;
    return Params();
  }


  template <typename T, typename Tag, typename Base>
  cgal_bgl_named_params<T,Tag,Base>
  inline no_parameters(cgal_bgl_named_params<T,Tag,Base>)
  {
    typedef cgal_bgl_named_params<T,Tag,Base> Params;
    return Params();
  }

// define free functions for named parameters
#define CGAL_add_named_parameter(X, Y, Z)        \
  template <typename K>                        \
  cgal_bgl_named_params<K, internal_np::X>                  \
  Z(K const& p)                                \
  {                                            \
    typedef cgal_bgl_named_params<K, internal_np::X> Params;\
    return Params(p);                          \
  }
#include <CGAL/boost/graph/parameters_interface.h>
#include <CGAL/boost/graph/boost_parameters_interface.h>
#undef CGAL_add_named_parameter

  } // namespace parameters

} //namespace CGAL

// partial specializations hate inheritance and we need to repeat
// those here. this is rather fragile.
namespace boost {

#if BOOST_VERSION < 105100
template <class Tag1, class Tag2, class T1, class Base>
inline
typename property_value< CGAL::cgal_bgl_named_params<T1,Tag1,Base>, Tag2>::type
get_param(const CGAL::cgal_bgl_named_params<T1,Tag1,Base>& p, Tag2 tag2)
{
  enum { match = detail::same_property<Tag1,Tag2>::value };
  typedef typename
    boost::property_value< CGAL::cgal_bgl_named_params<T1,Tag1,Base>, Tag2>::type T2;
  T2* t2 = 0;
  typedef detail::property_value_dispatch<match> Dispatcher;
  return Dispatcher::const_get_value(p, t2, tag2);
}
#endif


template <typename T, typename Tag, typename Base, typename Def>
struct lookup_named_param_def<Tag, CGAL::cgal_bgl_named_params<T, Tag, Base>, Def> {
  typedef T type;
  static const type& get(const bgl_named_params<T, Tag, Base>& p, const Def&) {
    return p.m_value;
  }
};

template <typename Tag1, typename T, typename Tag, typename Base, typename Def>
struct lookup_named_param_def<Tag1, CGAL::cgal_bgl_named_params<T, Tag, Base>, Def> {
  typedef typename lookup_named_param_def<Tag1, Base, Def>::type type;
  static const type& get(const bgl_named_params<T, Tag, Base>& p, const Def& def) {
    return lookup_named_param_def<Tag1, Base, Def>::get(p.m_base, def);
  }
};
} // boost

#include <CGAL/enable_warnings.h>

#endif // CGAL_BOOST_GRAPH_NAMED_FUNCTION_PARAMS_HPP
