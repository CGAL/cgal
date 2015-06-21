// Copyright (c) 2013 CNRS and LIRIS' Establishments (France).
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
// $URL:
// 
//
// Author(s)     : Pierre Talbot

#ifndef CGAL_BOOST_GRAPH_PROPERTIES_COMBINATORIAL_MAP_H
#define CGAL_BOOST_GRAPH_PROPERTIES_COMBINATORIAL_MAP_H

#include <CGAL/boost/graph/properties.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Linear_cell_complex.h>

#define CGAL_LCC_ARGS unsigned int d_, unsigned int ambient_dim, \
             class Traits_, \
             class Items_, \
             class Alloc_, \
             template<unsigned int, class, class, class> \
             class CMap

#define CGAL_LCC_TYPE CGAL::Linear_cell_complex<d_, ambient_dim, Traits_, Items_, Alloc_, CMap > 

namespace CGAL {

template<class LCC>
class LCC_edge_weight_map : public boost::put_get_helper<double, LCC>
{ 
public:

  typedef boost::readable_property_map_tag                                category;
  typedef double                                                          value_type;
  typedef double                                                          reference;
  typedef typename boost::graph_traits<LCC const>::edge_descriptor key_type;

  LCC_edge_weight_map(LCC const& ) {}

  reference operator[](key_type const& e) const
  {
    return CGAL::squared_distance(LCC::point(e), LCC::point(e->opposite()));
  }
};

template <class CMap>
class CMap_dart_index_map_external 
  : public boost::put_get_helper<std::size_t, CMap_dart_index_map_external<CMap> >
{
  
public:

  typedef boost::readable_property_map_tag                                category;
  typedef std::size_t                                                     value_type;
  typedef std::size_t                                                     reference;
  typedef typename boost::graph_traits<CMap const>::edge_descriptor       key_type;

  CMap_dart_index_map_external() : map()
  {}

  CMap_dart_index_map_external(CMap const& cm) 
    : map()
  {
    CMap &cmap = const_cast<CMap&>(cm);
    typedef typename CMap::template One_dart_per_cell_range<1>::iterator Iter;
    Iter b=cmap.template one_dart_per_cell<1>().begin();
    Iter e=cmap.template one_dart_per_cell<1>().end();
    for(value_type i=0; b != e; ++b, ++i)
    {
      map[b] = i;
      ++i;
      map[b->opposite()] = i;
    }
    typedef typename CMap::Dart_range::iterator DIter;
    DIter bi = cmap.darts().begin();
    DIter ei = cmap.darts().end();
    for(; bi != ei; ++bi)
    {
      value_type v1 = map[bi];
      value_type v2 = map[opposite_edge(bi, cm)];
      CGAL_assertion(v1 - v2 == 1 || v1 - v2 == -1);
    }
  }

  reference operator[](key_type const& e) const { return map[e]; }
  
private:
  
  CGAL::Unique_hash_map<key_type,std::size_t> map ;
};

template<class LCC>
class LCC_vertex_point_map : public boost::put_get_helper< typename LCC::Point&, LCC_vertex_point_map<LCC> >
{
public:

  typedef typename LCC::Point Point;
  
  typedef boost::lvalue_property_map_tag category;
  typedef Point value_type;
  typedef Point& reference;
  typedef typename boost::graph_traits<LCC>::vertex_descriptor key_type;

  LCC_vertex_point_map(LCC& ) {}

  reference operator[](key_type const& v) const { return LCC::point(v); }
};

template<class LCC>
class LCC_vertex_point_const_map : public boost::put_get_helper< typename LCC::Point const&, LCC_vertex_point_const_map<LCC> >
{
public:

  typedef typename LCC::Point Point;
  
  typedef boost::readable_property_map_tag category;
  typedef Point value_type;
  typedef Point const& reference;
  typedef typename boost::graph_traits<LCC const>::vertex_descriptor key_type;

  LCC_vertex_point_const_map(LCC const&) {}

  reference operator[](key_type const& v) const { return LCC::point(v); }
};


template <class CMap, typename Tag>
struct CMap_property_map
{};

template <class CMap>
struct CMap_property_map<CMap, vertex_external_index_t> 
{
  typedef CMap_dart_index_map_external<CMap> type;
  typedef CMap_dart_index_map_external<CMap> const_type;
};

template <class CMap>
struct CMap_property_map<CMap, edge_external_index_t> 
{
  typedef CMap_dart_index_map_external<CMap> type;
  typedef CMap_dart_index_map_external<CMap> const_type;
};

template <class CMap>
struct CMap_property_map<CMap, boost::edge_index_t>
{
  typedef CMap_dart_index_map_external<CMap> type;
  typedef CMap_dart_index_map_external<CMap> const_type;
};

template <class CMap>
struct CMap_property_map<CMap, boost::vertex_index_t>
{
  typedef CMap_dart_index_map_external<CMap> type;
  typedef CMap_dart_index_map_external<CMap> const_type;
};

// For LCC
template <class LCC>
struct CMap_property_map<LCC, boost::edge_weight_t>
{
  typedef LCC_edge_weight_map<LCC> type;
  typedef LCC_edge_weight_map<LCC> const_type;
};

template <class LCC>
struct CMap_property_map<LCC, vertex_point_t>
{
  typedef LCC_vertex_point_map<LCC> type;
  typedef LCC_vertex_point_const_map<LCC> const_type;
};

template <class Type, class Tag>
struct CMap_property_map_helper
{
  typedef typename CGAL::CMap_property_map<Type, Tag> map_gen;
  typedef typename map_gen::type       type;
  typedef typename map_gen::const_type const_type;
};

} // namespace CGAL

namespace boost{

// LCC.
template<CGAL_LCC_ARGS, class Tag>
struct property_map<CGAL_LCC_TYPE, Tag> : CGAL::CMap_property_map_helper<CGAL_LCC_TYPE, Tag>
{};

// This partial specialization shouldn't be needed but is due to a bug in Boost 1.51.
template<CGAL_LCC_ARGS, class Tag>
struct property_map<const CGAL_LCC_TYPE, Tag> : CGAL::CMap_property_map_helper<CGAL_LCC_TYPE, Tag>
{};

template<CGAL_LCC_ARGS>
CGAL::CMap_dart_index_map_external<CGAL_LCC_TYPE> get(CGAL::edge_external_index_t, CGAL_LCC_TYPE const& cmap) 
{
  return CGAL::CMap_dart_index_map_external<CGAL_LCC_TYPE>(cmap);
}

template<CGAL_LCC_ARGS>
CGAL::CMap_dart_index_map_external<CGAL_LCC_TYPE> get(CGAL::vertex_external_index_t, CGAL_LCC_TYPE const& cmap) 
{
  return CGAL::CMap_dart_index_map_external<CGAL_LCC_TYPE>(cmap);
}

// LCC get.
template <CGAL_LCC_ARGS>
CGAL::LCC_edge_weight_map<CGAL_LCC_TYPE> get(edge_weight_t, CGAL_LCC_TYPE const& lcc) 
{
  return CGAL::LCC_edge_weight_map<CGAL_LCC_TYPE>(lcc);
}

template <CGAL_LCC_ARGS>
CGAL::CMap_dart_index_map_external<CGAL_LCC_TYPE> get(edge_index_t, CGAL_LCC_TYPE const& lcc) 
{
  return CGAL::CMap_dart_index_map_external<CGAL_LCC_TYPE>(lcc);
}

template <CGAL_LCC_ARGS>
CGAL::CMap_dart_index_map_external<CGAL_LCC_TYPE> get(vertex_index_t, CGAL_LCC_TYPE const& lcc) 
{
  return CGAL::CMap_dart_index_map_external<CGAL_LCC_TYPE>(lcc);
}

template <CGAL_LCC_ARGS>
CGAL::LCC_vertex_point_map<CGAL_LCC_TYPE> get(CGAL::vertex_point_t, CGAL_LCC_TYPE& lcc) 
{
  return CGAL::LCC_vertex_point_map<CGAL_LCC_TYPE>(lcc);
}

template <CGAL_LCC_ARGS>
CGAL::LCC_vertex_point_const_map<CGAL_LCC_TYPE> get(CGAL::vertex_point_t, CGAL_LCC_TYPE const& lcc) 
{
  return CGAL::LCC_vertex_point_const_map<CGAL_LCC_TYPE>(lcc);
}

template<CGAL_LCC_ARGS, class PropertyTag, class Key>
typename property_traits<typename property_map<CGAL_LCC_TYPE, PropertyTag>::type>::reference
get(PropertyTag p, CGAL_LCC_TYPE& g, const Key& key) 
{
  return get(get(p, g), key);
}

template<CGAL_LCC_ARGS, class PropertyTag, class Key>
typename property_traits<typename property_map<CGAL_LCC_TYPE, PropertyTag>::const_type>::reference
get(PropertyTag p, CGAL_LCC_TYPE const& g, const Key& key) 
{
  return get(get(p, g), key);
}

template<CGAL_LCC_ARGS, class PropertyTag, class Key,class Value>
void put(PropertyTag p, CGAL_LCC_TYPE& g, const Key& key, const Value& value)
{
  typedef typename property_map<CGAL_LCC_TYPE, PropertyTag>::type Map;
  Map pmap = get(p, g);
  put(pmap, key, value);
}

} // namespace boost

#undef CGAL_LCC_ARGS
#undef CGAL_LCC_TYPE


#endif // CGAL_BOOST_GRAPH_PROPERTIES_COMBINATORIAL_MAP_H
