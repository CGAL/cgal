//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// This file is part of the Boost Graph Library
//
// You should have received a copy of the License Agreement for the
// Boost Graph Library along with the software; see the file LICENSE.
// If not, contact Office of Research, University of Notre Dame, Notre
// Dame, IN 46556.
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.

//=======================================================================
// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri, Fernando Cacciola



#ifndef CGAL_BOOST_GRAPH_NAMED_FUNCTION_PARAMS_H
#define CGAL_BOOST_GRAPH_NAMED_FUNCTION_PARAMS_H

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


namespace CGAL {
#define CGAL_add_pmp_parameter(X, Y, Z)            \
  enum X { Y };                                    \

#include <CGAL/boost/graph/parameters_interface.h>
#undef CGAL_add_pmp_parameter

  template <typename T, typename Tag, typename Base = boost::no_property>
  struct cgal_bgl_named_params : boost::bgl_named_params<T, Tag, Base>
  {
    typedef boost::bgl_named_params<T, Tag, Base> base;
    typedef cgal_bgl_named_params self;

    cgal_bgl_named_params(T v = T()) : base(v) {}
    cgal_bgl_named_params(T v, const Base& b) : base(v, b) {}


#define CGAL_add_pmp_parameter(X, Y, Z)              \
  template<typename K>                               \
  cgal_bgl_named_params<K, X, self>                  \
  Z(const K& k) const                                \
  {                                                  \
    typedef cgal_bgl_named_params<K, X, self> Params;\
    return Params(k, *this);                         \
  }
#include <CGAL/boost/graph/parameters_interface.h>
#undef CGAL_add_pmp_parameter

#define CGAL_add_pmp_parameter(X, Y, Z)                      \
  template<typename K>                                       \
  cgal_bgl_named_params<K, boost::X, self>                   \
  Z(const K& k) const                                        \
  {                                                          \
    typedef cgal_bgl_named_params<K, boost::X, self> Params; \
    return Params(k, *this);                                 \
  }
#include <CGAL/boost/graph/boost_parameters_interface.h>
#undef CGAL_add_pmp_parameter
    template <typename IndexMap>
    cgal_bgl_named_params<IndexMap, boost::edge_index_t, self>
    edge_index_map(const IndexMap& p) const 
    {
      typedef cgal_bgl_named_params<IndexMap, boost::edge_index_t, self> Params;
      return Params(p, *this);
    }

      template <typename IndexMap>
    cgal_bgl_named_params<IndexMap, boost::halfedge_index_t, self>
    halfedge_index_map(const IndexMap& p) const 
    {
      typedef cgal_bgl_named_params<IndexMap, boost::halfedge_index_t, self> Params;
      return Params(p, *this);
    }

    template <typename Visitor>
    cgal_bgl_named_params<Visitor, boost::graph_visitor_t, self>
    visitor(const Visitor& p) const 
    {
      typedef cgal_bgl_named_params<Visitor, boost::graph_visitor_t, self> Params;
      return Params(p, *this);
    }
  };

  namespace parameters {
  
  template <typename IndexMap>
  cgal_bgl_named_params<IndexMap, boost::halfedge_index_t>
  halfedge_index_map(IndexMap const& p) 
  {
    typedef cgal_bgl_named_params<IndexMap, boost::halfedge_index_t> Params;
    return Params(p);
  }
#define CGAL_add_pmp_parameter(X, Y, Z)               \
  template <typename K>                               \
  cgal_bgl_named_params<K, boost::X>                  \
  Z(K const& p)                                       \
  {                                                   \
    typedef cgal_bgl_named_params<K, boost::X> Params;\
    return Params(p);                                 \
  }
#include <CGAL/boost/graph/boost_parameters_interface.h>
#undef CGAL_add_pmp_parameter

#define CGAL_add_pmp_parameter(X, Y, Z)        \
  template <typename K>                        \
  cgal_bgl_named_params<K, X>                  \
  Z(K const& p)                                \
  {                                            \
    typedef cgal_bgl_named_params<K, X> Params;\
    return Params(p);                          \
  }
#include <CGAL/boost/graph/parameters_interface.h>
#undef CGAL_add_pmp_parameter
  template <typename IndexMap>
  cgal_bgl_named_params<IndexMap, boost::edge_index_t>
  edge_index_map(IndexMap const& pmap) 
  {
    typedef cgal_bgl_named_params<IndexMap, boost::edge_index_t> Params;
    return Params(pmap);
  }

  template <typename Visitor>
  cgal_bgl_named_params<Visitor, boost::graph_visitor_t>
  visitor(const Visitor& p) 
  {
    typedef cgal_bgl_named_params<Visitor, boost::graph_visitor_t> Params;
    return Params(p);
  }
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

#endif // CGAL_BOOST_GRAPH_NAMED_FUNCTION_PARAMS_HPP
