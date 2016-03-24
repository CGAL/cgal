// Copyright (c) 2016  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Guillaume Damiand

#ifndef CGAL_BOOST_GRAPH_PROPERTIES_LINEAR_CELL_COMPLEX_H
#define CGAL_BOOST_GRAPH_PROPERTIES_LINEAR_CELL_COMPLEX_H

#include <CGAL/boost/graph/properties.h>
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Unique_hash_map.h>

#define CGAL_LCC_ARGS unsigned int d_, unsigned int ambient_dim,        \
             class Traits_, \
             class Items_, \
             class Alloc_, \
             template<unsigned int, class, class, class, class> \
             class CMap, \
             class Storage_

#define CGAL_LCC_TYPE CGAL::Linear_cell_complex<d_, ambient_dim, Traits_, Items_, Alloc_, CMap , Storage_>

namespace CGAL {

template<class LCC>
class LCC_edge_weight_map : public boost::put_get_helper<double, LCC_edge_weight_map<LCC> >
{ 
public:

  typedef boost::readable_property_map_tag                   category;
  typedef double                                             value_type;
  typedef double                                             reference;
  typedef typename boost::graph_traits<LCC>::edge_descriptor key_type;

  LCC_edge_weight_map(LCC const& lcc) : m_lcc(lcc) {}

  //  reference operator[](key_type const& e)
  //{ return CGAL::squared_distance(m_lcc.point(e), m_lcc.point(m_lcc.other_extremity(e))); }
  reference operator[](key_type const& e) const
  { return CGAL::squared_distance(m_lcc.point(e), m_lcc.point(m_lcc.other_extremity(e))); }

protected:
  LCC const& m_lcc;
};

template <class CMap>
class CMap_dart_index_map_external :
    public boost::put_get_helper<std::size_t&, CMap_dart_index_map_external<CMap> >
{
  
public:

  typedef boost::read_write_property_map_tag                  category;
  typedef std::size_t                                         value_type;
  typedef std::size_t&                                        reference;
  typedef typename boost::graph_traits<CMap>::edge_descriptor key_type;

  CMap_dart_index_map_external()
  {}

  CMap_dart_index_map_external(CMap const& cm) 
  {
    CMap &cmap = const_cast<CMap&>(cm);
    typedef typename CMap::template One_dart_per_cell_range<1>::iterator Iter;
    Iter b=cmap.template one_dart_per_cell<1>().begin();
    Iter e=cmap.template one_dart_per_cell<1>().end();
    for(value_type i=0; b != e; ++b, ++i)
    {
      map[b] = i;
      ++i;
      map[cmap.template beta<2>(b)] = i;
    }
    typedef typename CMap::Dart_range::iterator DIter;
    DIter bi = cmap.darts().begin();
    DIter ei = cmap.darts().end();
    for(; bi != ei; ++bi)
    {
      value_type v1 = map[bi];
      value_type v2 = map[cmap.template beta<2>(bi)];
      CGAL_assertion(v1 - v2 == 1 || v1 - v2 == -1);
    }
  }

  reference operator[](key_type const& e) const
  { return const_cast<CGAL::Unique_hash_map<key_type,std::size_t>&>(map)[e]; }
  
private:  
  CGAL::Unique_hash_map<key_type,std::size_t> map ;
};

template <class CMap>
class CMap_vertex_index :
    public boost::put_get_helper<std::size_t&, CMap_vertex_index<CMap> >
{
public:
  typedef boost::read_write_property_map_tag                    category;
  typedef std::size_t                                           value_type;
  typedef std::size_t&                                          reference;
  typedef typename boost::graph_traits<CMap>::vertex_descriptor key_type;

  CMap_vertex_index()
  {}

  CMap_vertex_index(CMap const& cm) 
  {
    CMap &cmap = const_cast<CMap&>(cm);
    typedef typename CMap::template Attribute_range<0>::type::iterator Iter;
    Iter b=cmap.template attributes<0>().begin();
    Iter e=cmap.template attributes<0>().end();
    for(value_type i=0; b != e; ++b, ++i)
    {
      map[b] = i;
    }
  }

  reference operator[](key_type const& e) const
  { return const_cast<CGAL::Unique_hash_map<key_type,std::size_t>&>(map)[e]; }
  
private:  
  CGAL::Unique_hash_map<key_type,std::size_t> map ;
};

template <class CMap>
class CMap_face_index :
    public boost::put_get_helper<std::size_t&, CMap_face_index<CMap> >
{
public:
  typedef boost::read_write_property_map_tag                    category;
  typedef std::size_t                                           value_type;
  typedef std::size_t&                                          reference;
  typedef typename boost::graph_traits<CMap>::face_descriptor key_type;

  CMap_face_index()
  {}

  CMap_face_index(CMap const& cm)
  {
    CMap &cmap = const_cast<CMap&>(cm);
    typedef typename CMap::template Attribute_range<2>::type::iterator Iter;
    Iter b=cmap.template attributes<2>().begin();
    Iter e=cmap.template attributes<2>().end();
    for(value_type i=0; b != e; ++b, ++i)
    {
      map[b] = i;
    }
  }

  reference operator[](key_type const& e) const
  { return const_cast<CGAL::Unique_hash_map<key_type,std::size_t>&>(map)[e]; }

private:
  CGAL::Unique_hash_map<key_type,std::size_t> map ;
};


template<class LCC>
class LCC_vertex_point_map :
    public boost::put_get_helper< typename LCC::Point&, LCC_vertex_point_map<LCC> >
{
public:

  typedef typename LCC::Point Point;
  
  typedef boost::lvalue_property_map_tag category;
  typedef Point value_type;
  typedef Point& reference;
  typedef typename boost::graph_traits<LCC>::vertex_descriptor key_type;

  LCC_vertex_point_map(LCC& lcc) : m_lcc(lcc) {}

  reference operator[](key_type const& v) const
  { return m_lcc.point_of_vertex_attribute(v); }

protected:
  LCC & m_lcc;
};

template<class LCC>
class LCC_vertex_point_const_map :
    public boost::put_get_helper< typename LCC::Point const&, LCC_vertex_point_const_map<LCC> >
{
public:
  typedef typename LCC::Point Point;
  
  typedef boost::readable_property_map_tag category;
  typedef Point value_type;
  typedef Point const& reference;
  typedef typename boost::graph_traits<LCC const>::vertex_descriptor key_type;

  LCC_vertex_point_const_map(LCC const& lcc) : m_lcc(lcc) {}

  reference operator[](key_type const& v) const
  { return m_lcc.point_of_vertex_attribute(v); }

protected:
  const LCC & m_lcc;
};

template <class CMap, typename Tag>
struct CMap_property_map
{};

template <class CMap>
struct CMap_property_map<CMap, boost::vertex_index_t>
{
  typedef CMap_vertex_index<CMap> type;
  typedef CMap_vertex_index<CMap> const_type;
};

template <class CMap>
struct CMap_property_map<CMap, boost::vertex_external_index_t> 
{
  typedef CMap_vertex_index<CMap> type;
  typedef CMap_vertex_index<CMap> const_type;
};

template <class CMap>
struct CMap_property_map<CMap, boost::edge_index_t>
{
  typedef CMap_dart_index_map_external<CMap> type;
  typedef CMap_dart_index_map_external<CMap> const_type;
};

template <class CMap>
struct CMap_property_map<CMap, boost::edge_external_index_t> 
{
  typedef CMap_dart_index_map_external<CMap> type;
  typedef CMap_dart_index_map_external<CMap> const_type;
};

template <class CMap>
struct CMap_property_map<CMap, boost::face_index_t>
{
  typedef CMap_face_index<CMap> type;
  typedef CMap_face_index<CMap> const_type;
};

template <class CMap>
struct CMap_property_map<CMap, boost::face_external_index_t> 
{
  typedef CMap_face_index<CMap> type;
  typedef CMap_face_index<CMap> const_type;
};

template <class CMap>
struct CMap_property_map<CMap, boost::halfedge_index_t>
{
  typedef CMap_dart_index_map_external<CMap> type;
  typedef CMap_dart_index_map_external<CMap> const_type;
};

template <class CMap>
struct CMap_property_map<CMap, boost::halfedge_external_index_t>
{
  typedef CMap_dart_index_map_external<CMap> type;
  typedef CMap_dart_index_map_external<CMap> const_type;
};

template <class LCC>
struct CMap_property_map<LCC, boost::edge_weight_t>
{
  typedef LCC_edge_weight_map<LCC> type;
  typedef LCC_edge_weight_map<LCC> const_type;
};

template <class LCC>
struct CMap_property_map<LCC, boost::vertex_point_t>
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
struct edge_property_type<CGAL_LCC_TYPE>
{
  typedef edge_weight_t type;
};

template<CGAL_LCC_ARGS>
struct vertex_property_type<CGAL_LCC_TYPE>
{
  typedef CGAL::vertex_point_t type;
};

template<CGAL_LCC_ARGS>
struct vertex_property_type<const CGAL_LCC_TYPE>
{
  typedef CGAL::vertex_point_t type;
};
} // namespace boost

namespace CGAL{
  
template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::halfedge_index_t >::const_type
get(boost::halfedge_index_t, CGAL_LCC_TYPE const& cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::halfedge_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::halfedge_external_index_t >::const_type
get(boost::halfedge_external_index_t, CGAL_LCC_TYPE const& cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::halfedge_external_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_index_t >::const_type
get(boost::vertex_index_t, CGAL_LCC_TYPE const& cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::vertex_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_external_index_t >::const_type
get(boost::vertex_external_index_t, CGAL_LCC_TYPE const& cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::vertex_external_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::edge_index_t >::const_type
get(boost::edge_index_t, CGAL_LCC_TYPE const& cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::edge_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::edge_external_index_t >::const_type
get(boost::edge_external_index_t, CGAL_LCC_TYPE const& cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::edge_external_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::face_index_t >::const_type
get(boost::face_index_t, CGAL_LCC_TYPE const& cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::face_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}
  
template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::face_external_index_t >::const_type
get(boost::face_external_index_t, CGAL_LCC_TYPE const& cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::face_external_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_point_t >::const_type
get(boost::vertex_point_t, CGAL_LCC_TYPE const& cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::vertex_point_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::edge_weight_t >::const_type
get(boost::edge_weight_t, CGAL_LCC_TYPE const& cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::edge_weight_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}
 
// the same blurb for non-const
template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::halfedge_index_t >::type
get(boost::halfedge_index_t, CGAL_LCC_TYPE & cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::halfedge_index_t>::
    type(cmap);
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::halfedge_external_index_t >::type
get(boost::halfedge_external_index_t, CGAL_LCC_TYPE & cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::halfedge_external_index_t>::
    type(cmap);
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_index_t >::type
get(boost::vertex_index_t, CGAL_LCC_TYPE & cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::vertex_index_t>::
    type(cmap);
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_external_index_t >::type
get(boost::vertex_external_index_t, CGAL_LCC_TYPE & cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::vertex_external_index_t>::
    type(cmap);
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::edge_index_t >::type
get(boost::edge_index_t, CGAL_LCC_TYPE & cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::edge_index_t>::
    type(cmap);
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::edge_external_index_t >::type
get(boost::edge_external_index_t, CGAL_LCC_TYPE & cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::edge_external_index_t>::
    type(cmap);
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::face_external_index_t >::type
get(boost::face_external_index_t, CGAL_LCC_TYPE & cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::face_external_index_t>::
    type(cmap);
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::face_index_t >::type
get(boost::face_index_t, CGAL_LCC_TYPE & cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::face_index_t>::
    type(cmap);
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_point_t >::type
get(boost::vertex_point_t, CGAL_LCC_TYPE & cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::vertex_point_t>::
    type(cmap);
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::edge_weight_t >::type
get(boost::edge_weight_t, CGAL_LCC_TYPE & cmap) 
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::edge_weight_t>::
    type(cmap);
}

///// 
template<CGAL_LCC_ARGS, class PropertyTag, class Key>
typename boost::property_traits<typename boost::property_map<CGAL_LCC_TYPE, PropertyTag>::type>::reference
get(PropertyTag p, CGAL_LCC_TYPE& g, const Key& key) 
{
  return get(get(p, g), key);
}

template<CGAL_LCC_ARGS, class PropertyTag, class Key>
typename boost::property_traits<typename boost::property_map<CGAL_LCC_TYPE, PropertyTag>::const_type>::reference
get(PropertyTag p, CGAL_LCC_TYPE const& g, const Key& key) 
{
  return get(get(p, g), key);
}

template<CGAL_LCC_ARGS, class PropertyTag, class Key,class Value>
void put(PropertyTag p, CGAL_LCC_TYPE& g, const Key& key, const Value& value)
{
  typedef typename boost::property_map<CGAL_LCC_TYPE, PropertyTag>::type Map;
  Map pmap = get(p, g);
  put(pmap, key, value);
}

} // namespace CGAL

#undef CGAL_LCC_ARGS
#undef CGAL_LCC_TYPE

#endif // CGAL_BOOST_GRAPH_PROPERTIES_LINEAR_CELL_COMPLEX_H

