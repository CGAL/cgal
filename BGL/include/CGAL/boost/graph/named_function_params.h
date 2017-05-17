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

#define CGAL_BGL_NP_TEMPLATE_PARAMETERS T, typename Tag, typename Base
#define CGAL_BGL_NP_CLASS CGAL::cgal_bgl_named_params<T,Tag,Base>

namespace CGAL {
namespace internal_np{

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

  //helper classes
  template<typename PolygonMesh, typename PropertyTag>
  class property_map_selector
  {
  public:
    typedef typename boost::graph_has_property<PolygonMesh, PropertyTag>::type Has_internal_pmap;
    typedef typename boost::mpl::if_c< Has_internal_pmap::value
                            , typename boost::property_map<PolygonMesh, PropertyTag>::type
                            , typename boost::cgal_no_property::type
    >::type type;
    typedef typename boost::mpl::if_c< Has_internal_pmap::value
                            , typename boost::property_map<PolygonMesh, PropertyTag>::const_type
                            , typename boost::cgal_no_property::const_type
    >::type const_type;

    type get_pmap(const PropertyTag& p, PolygonMesh& pmesh)
    {
      return get_impl(p, pmesh, Has_internal_pmap());
    }

    const_type get_const_pmap(const PropertyTag& p, const PolygonMesh& pmesh)
    {
      return get_const_pmap_impl(p, pmesh, Has_internal_pmap());
    }

  private:
    type get_impl(const PropertyTag&, PolygonMesh&, CGAL::Tag_false)
    {
      return type(); //boost::cgal_no_property::type
    }
    type get_impl(const PropertyTag& p, PolygonMesh& pmesh, CGAL::Tag_true)
    {
      return get(p, pmesh);
    }

    const_type get_const_pmap_impl(const PropertyTag&
                                 , const PolygonMesh&, CGAL::Tag_false)
    {
      return const_type(); //boost::cgal_no_property::type
    }
    const_type get_const_pmap_impl(const PropertyTag& p
                                 , const PolygonMesh& pmesh, CGAL::Tag_true)
    {
      return get(p, pmesh);
    }
  };

  template<typename PolygonMesh, typename PropertyTag>
  typename property_map_selector<PolygonMesh, PropertyTag>::type
  get_property_map(const PropertyTag& p, PolygonMesh& pmesh)
  {
    property_map_selector<PolygonMesh, PropertyTag> pms;
    return pms.get_pmap(p, pmesh);
  }

  template<typename PolygonMesh, typename PropertyTag>
  typename property_map_selector<PolygonMesh, PropertyTag>::const_type
  get_const_property_map(const PropertyTag& p, const PolygonMesh& pmesh)
  {
    property_map_selector<PolygonMesh, PropertyTag> pms;
    return pms.get_const_pmap(p, pmesh);
  }
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
