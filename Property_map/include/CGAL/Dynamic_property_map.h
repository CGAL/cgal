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
// Author(s)     : Andreas Fabri

#ifndef CGAL_DYNAMIC_PROPERTY_MAP_H
#define CGAL_DYNAMIC_PROPERTY_MAP_H

#include <CGAL/boost/graph/properties.h>
#include <CGAL/property_map.h>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

namespace CGAL {

namespace internal {

template <typename K, typename V>
struct Dynamic_property_map {

  typedef K key_type;
  typedef V value_type;
  typedef value_type& reference;
  typedef boost::read_write_property_map_tag  category;

  Dynamic_property_map(const V& default_value = V())
    : map_(new Map()), default_value_(default_value)
  {}

  void clear()
  {
    map_ = boost::shared_ptr<Map>(0);
  }


  friend value_type get(const Dynamic_property_map& m, const key_type& k)
  {
    typename Map::const_iterator it = m.map_->find(k);
    if(it == m.map_->end()){
      return m.default_value();
    }
    return it->second;
  }


  friend void put(Dynamic_property_map& m, const key_type& k, const value_type& v)
  {
    if(v != m.default_value()){
      (*(m.map_))[k] = v;
    }
  }

  
  const V& default_value() const
  {
    return default_value_;
  }


  typedef boost::unordered_map<K,V> Map;
  boost::shared_ptr<Map> map_;
  V default_value_;
};
} // namespace internal

} // namespace CGAL

namespace boost {

template < typename T >
struct vertex_property_t;

template <typename G, typename T>
struct property_map<G, boost::vertex_property_t<T> >
{
  typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
  typedef CGAL::internal::Dynamic_property_map<vertex_descriptor,T> type;
  typedef type const_type;
};

template <class G, typename T>
struct property_map<G, boost::halfedge_property_t<T> >
{
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<halfedge_descriptor,T> type;
  typedef type const_type;
};

template <class G, typename T>
struct property_map<G, boost::edge_property_t<T> >
{
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<edge_descriptor,T> type;
  typedef type const_type;
};

template <class G, typename T>
struct property_map<G, boost::face_property_t<T> >
{
  typedef typename boost::graph_traits<G>::face_descriptor face_descriptor;
  typedef CGAL::internal::Dynamic_property_map<face_descriptor,T> type;
  typedef type const_type;
};
}

#endif // CGAL_DYNAMIC_PROPERTY_MAP_H
