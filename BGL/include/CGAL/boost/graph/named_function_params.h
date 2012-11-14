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
#include <boost/graph/properties.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/version.hpp>

// An explanation about the version hackery below: There is no real
// API to introduce custom properties to the Graph API and the
// internals have changed with Boost Version 1.51 and changes aren't
// backward compatible. To work around that we carry around two
// versions of cgal_bgl_named_params. One imitates the pre 1.51
// bgl_named_params, the newer one hooks into the API through
// inheritance and addition of the some partial specializations.
namespace CGAL {

  enum vertex_is_fixed_t             { vertex_is_fixed             } ;  
  enum set_cache_policy_t            { set_cache_policy            } ;
  enum get_cost_policy_t             { get_cost_policy             } ;
  enum get_cost_policy_params_t      { get_cost_policy_params      } ;
  enum get_placement_policy_t        { get_placement_policy        } ;
  enum get_placement_policy_params_t { get_placement_policy_params } ;
  
#if BOOST_VERSION >= 105100
  template <typename T, typename Tag, typename Base = boost::no_property>
  struct cgal_bgl_named_params : boost::bgl_named_params<T, Tag, Base>
  {
    typedef boost::bgl_named_params<T, Tag, Base> base;
    typedef cgal_bgl_named_params self;

    cgal_bgl_named_params(T v = T()) : base(v) {}
    cgal_bgl_named_params(T v, const Base& b) : base(v, b) {}

    template <typename IndexMap>
    cgal_bgl_named_params<IndexMap, boost::vertex_index_t, self>
    vertex_index_map(const IndexMap& p) const 
    {
      typedef cgal_bgl_named_params<IndexMap, boost::vertex_index_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename PointMap>
    cgal_bgl_named_params<PointMap, vertex_point_t, self>
    vertex_point_map(const PointMap& p) const 
    {
      typedef cgal_bgl_named_params<PointMap, vertex_point_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename IsFixedMap>
    cgal_bgl_named_params<IsFixedMap, vertex_is_fixed_t, self>
    vertex_is_fixed_map(const IsFixedMap& p) const 
    {
      typedef cgal_bgl_named_params<IsFixedMap, vertex_is_fixed_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename IndexMap>
    cgal_bgl_named_params<IndexMap, boost::edge_index_t, self>
    edge_index_map(const IndexMap& p) const 
    {
      typedef cgal_bgl_named_params<IndexMap, boost::edge_index_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename IsBorderMap>
    cgal_bgl_named_params<IsBorderMap, edge_is_border_t, self>
    edge_is_border_map(const IsBorderMap& p) const 
    {
      typedef cgal_bgl_named_params<IsBorderMap, edge_is_border_t, self> Params;
      return Params(p, *this);
    }

    template <typename Visitor>
    cgal_bgl_named_params<Visitor, boost::graph_visitor_t, self>
    visitor(const Visitor& p) const 
    {
      typedef cgal_bgl_named_params<Visitor, boost::graph_visitor_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename SetCache>
    cgal_bgl_named_params<SetCache, set_cache_policy_t, self>
    set_cache(const SetCache& p) const 
    {
      typedef cgal_bgl_named_params<SetCache, set_cache_policy_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename GetCost>
    cgal_bgl_named_params<GetCost, get_cost_policy_t, self>
    get_cost(const GetCost& p) const 
    {
      typedef cgal_bgl_named_params<GetCost, get_cost_policy_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename GetCostParams>
    cgal_bgl_named_params<GetCostParams, get_cost_policy_params_t, self>
    get_cost_params(const GetCostParams& p) const 
    {
      typedef cgal_bgl_named_params<GetCostParams, get_cost_policy_params_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename GetPlacement>
    cgal_bgl_named_params<GetPlacement, get_placement_policy_t, self>
    get_placement(const GetPlacement& p) const 
    {
      typedef cgal_bgl_named_params<GetPlacement, get_placement_policy_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename GetPlacementParams>
    cgal_bgl_named_params<GetPlacementParams, get_placement_policy_params_t, self>
    get_placement_params(const GetPlacementParams& p) const 
    {
      typedef cgal_bgl_named_params<GetPlacementParams, get_placement_policy_params_t, self> Params;
      return Params(p, *this);
    }
  };
#else
  template <typename T, typename Tag, typename Base = boost::no_property>
  struct cgal_bgl_named_params : public Base
  {
    typedef cgal_bgl_named_params self;
    typedef Base next_type;
    typedef Tag tag_type;
    typedef T value_type;

    cgal_bgl_named_params() { }

    cgal_bgl_named_params(T v) : m_value(v) { }
    cgal_bgl_named_params(T v, const Base& b) : Base(b), m_value(v) { }

    T m_value;

    template <typename IndexMap>
    cgal_bgl_named_params<IndexMap, boost::vertex_index_t, self>
    vertex_index_map(const IndexMap& p) const 
    {
      typedef cgal_bgl_named_params<IndexMap, boost::vertex_index_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename PointMap>
    cgal_bgl_named_params<PointMap, vertex_point_t, self>
    vertex_point_map(const PointMap& p) const 
    {
      typedef cgal_bgl_named_params<PointMap, vertex_point_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename IsFixedMap>
    cgal_bgl_named_params<IsFixedMap, vertex_is_fixed_t, self>
    vertex_is_fixed_map(const IsFixedMap& p) const 
    {
      typedef cgal_bgl_named_params<IsFixedMap, vertex_is_fixed_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename IndexMap>
    cgal_bgl_named_params<IndexMap, boost::edge_index_t, self>
    edge_index_map(const IndexMap& p) const 
    {
      typedef cgal_bgl_named_params<IndexMap, boost::edge_index_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename IsBorderMap>
    cgal_bgl_named_params<IsBorderMap, edge_is_border_t, self>
    edge_is_border_map(const IsBorderMap& p) const 
    {
      typedef cgal_bgl_named_params<IsBorderMap, edge_is_border_t, self> Params;
      return Params(p, *this);
    }

    template <typename Visitor>
    cgal_bgl_named_params<Visitor, boost::graph_visitor_t, self>
    visitor(const Visitor& p) const 
    {
      typedef cgal_bgl_named_params<Visitor, boost::graph_visitor_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename SetCache>
    cgal_bgl_named_params<SetCache, set_cache_policy_t, self>
    set_cache(const SetCache& p) const 
    {
      typedef cgal_bgl_named_params<SetCache, set_cache_policy_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename GetCost>
    cgal_bgl_named_params<GetCost, get_cost_policy_t, self>
    get_cost(const GetCost& p) const 
    {
      typedef cgal_bgl_named_params<GetCost, get_cost_policy_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename GetCostParams>
    cgal_bgl_named_params<GetCostParams, get_cost_policy_params_t, self>
    get_cost_params(const GetCostParams& p) const 
    {
      typedef cgal_bgl_named_params<GetCostParams, get_cost_policy_params_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename GetPlacement>
    cgal_bgl_named_params<GetPlacement, get_placement_policy_t, self>
    get_placement(const GetPlacement& p) const 
    {
      typedef cgal_bgl_named_params<GetPlacement, get_placement_policy_t, self> Params;
      return Params(p, *this);
    }
    
    template <typename GetPlacementParams>
    cgal_bgl_named_params<GetPlacementParams, get_placement_policy_params_t, self>
    get_placement_params(const GetPlacementParams& p) const 
    {
      typedef cgal_bgl_named_params<GetPlacementParams, get_placement_policy_params_t, self> Params;
      return Params(p, *this);
    }
  };

  template <class Tag1, class Tag2, class T1, class Base>
  inline
  typename boost::property_value< cgal_bgl_named_params<T1,Tag1,Base>, Tag2>::type
  get_param(const cgal_bgl_named_params<T1,Tag1,Base>& p, Tag2 tag2)
  {
    enum { match = boost::detail::same_property<Tag1,Tag2>::value };
    typedef typename
      boost::property_value< cgal_bgl_named_params<T1,Tag1,Base>, Tag2>::type T2;
    T2* t2 = 0;
    typedef boost::detail::property_value_dispatch<match> Dispatcher;
    return Dispatcher::const_get_value(p, t2, tag2);
  }
#endif

  template <typename IndexMap>
  cgal_bgl_named_params<IndexMap, boost::vertex_index_t>
  vertex_index_map(IndexMap const& p) 
  {
    typedef cgal_bgl_named_params<IndexMap, boost::vertex_index_t> Params;
    return Params(p);
  }
  
  template <typename PointMap>
  cgal_bgl_named_params<PointMap, vertex_point_t>
  vertex_point_map(PointMap const& p) 
  {
    typedef cgal_bgl_named_params<PointMap, vertex_point_t> Params;
    return Params(p);
  }

  template <typename IsFixedMap>
  cgal_bgl_named_params<IsFixedMap, vertex_is_fixed_t>
  vertex_is_fixed_map(IsFixedMap const& p) 
  {
    typedef cgal_bgl_named_params<IsFixedMap, vertex_is_fixed_t> Params;
    return Params(p);
  }
  
  template <typename IndexMap>
  cgal_bgl_named_params<IndexMap, boost::edge_index_t>
  edge_index_map(IndexMap const& pmap) 
  {
    typedef cgal_bgl_named_params<IndexMap, boost::edge_index_t> Params;
    return Params(pmap);
  }
  
  template <typename IsBorderMap>
  cgal_bgl_named_params<IsBorderMap, edge_is_border_t>
  edge_is_border_map(IsBorderMap const& p) 
  {
    typedef cgal_bgl_named_params<IsBorderMap, edge_is_border_t> Params;
    return Params(p);
  }
  
  template <typename Visitor>
  cgal_bgl_named_params<Visitor, boost::graph_visitor_t>
  visitor(const Visitor& p) 
  {
    typedef cgal_bgl_named_params<Visitor, boost::graph_visitor_t> Params;
    return Params(p);
  }

  template <typename SetCache>
  cgal_bgl_named_params<SetCache, set_cache_policy_t>
  set_cache(const SetCache& p) 
  {
    typedef cgal_bgl_named_params<SetCache, set_cache_policy_t> Params;
    return Params(p);
  }
  
  template <typename GetCost>
  cgal_bgl_named_params<GetCost, get_cost_policy_t>
  get_cost(const GetCost& p)
  {
    typedef cgal_bgl_named_params<GetCost, get_cost_policy_t> Params;
     return Params(p);
  }

  template <typename GetCostParams>
  cgal_bgl_named_params<GetCostParams, get_cost_policy_params_t>
  get_cost_params(const GetCostParams& p)
  {
    typedef cgal_bgl_named_params<GetCostParams, get_cost_policy_params_t> Params;
    return Params(p);
  }
  
  template <typename GetPlacement>
  cgal_bgl_named_params<GetPlacement, get_placement_policy_t>
  get_placement(const GetPlacement& p) 
  {
    typedef cgal_bgl_named_params<GetPlacement, get_placement_policy_t> Params;
    return Params(p);
  }
  
  template <typename GetPlacementParams>
  cgal_bgl_named_params<GetPlacementParams, get_placement_policy_params_t>
  get_placement_params(const GetPlacementParams& p)
  {
    typedef cgal_bgl_named_params<GetPlacementParams, get_placement_policy_params_t> Params;
    return Params(p);
  }
} //namespace CGAL

#if BOOST_VERSION >= 105100
// partial specializations hate inheritance and we need to repeat
// those here. this is rather fragile.
namespace boost {
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
#endif

#endif // CGAL_BOOST_GRAPH_NAMED_FUNCTION_PARAMS_HPP
