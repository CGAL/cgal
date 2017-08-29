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

#include <boost/graph/graph_traits.hpp>
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
  typedef const value_type& reference;
  typedef boost::read_write_property_map_tag  category;

  Dynamic_property_map(const V& default_value = V())
    : map_(new Map()), default_value_(default_value)
  {}

  void clear()
  {
    if(map_){
      map_->clear();
    }
  }


  friend reference get(const Dynamic_property_map& m, const key_type& k)
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

template <typename T>
struct vertex_property_t
{
  vertex_property_t(const std::string s, const T& t = T())
    : s(s), t(t)
  {}
  std::string s;
  T t;
};


template <typename T>
struct halfedge_property_t
{
  halfedge_property_t(const std::string s, const T& t = T())
    : s(s), t(t)
  {}
  std::string s;
  T t;
};

template <typename T>
struct edge_property_t
{
  edge_property_t(const std::string s, const T& t = T())
    : s(s), t(t)
  {}
  std::string s;
  T t;
};


template <typename T>
struct face_property_t
{
  face_property_t(const std::string s, const T& t = T())
    : s(s), t(t)
  {}
  std::string s;
  T t;
};

template <typename G, typename Tag>
struct dynamic_property_map{};

template <typename G, typename T>
struct dynamic_property_map<G,vertex_property_t<T> >
{
  typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
  typedef CGAL::internal::Dynamic_property_map<vertex_descriptor,T> type;
  typedef type const_type;
};

template <typename G, typename T>
struct dynamic_property_map<G,halfedge_property_t<T> >
{
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<halfedge_descriptor,T> type;
  typedef type const_type;
};


template <typename G, typename T>
struct dynamic_property_map<G,edge_property_t<T> >
{
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<edge_descriptor,T> type;
  typedef type const_type;
};

template <typename G, typename T>
struct dynamic_property_map<G,face_property_t<T> >
{
  typedef typename boost::graph_traits<G>::face_descriptor face_descriptor;
  typedef CGAL::internal::Dynamic_property_map<face_descriptor,T> type;
  typedef type const_type;
};



template <typename T, typename G>
typename dynamic_property_map<G,vertex_property_t<T> >::const_type
add_property(vertex_property_t<T> prop, const G&)
{
  typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
  return internal::Dynamic_property_map<vertex_descriptor,T>(prop.t);
}

template <typename T, typename G>
typename dynamic_property_map<G,halfedge_property_t<T> >::const_type
add_property(halfedge_property_t<T> prop, const G&)
{
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  return internal::Dynamic_property_map<halfedge_descriptor,T>(prop.t);
}

template <typename T, typename G>
typename dynamic_property_map<G,edge_property_t<T> >::const_type
add_property(edge_property_t<T> prop, const G&)
{
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  return internal::Dynamic_property_map<edge_descriptor,T>(prop.t);
}

template <typename T, typename G>
typename dynamic_property_map<G,face_property_t<T> >::const_type
add_property(face_property_t<T> prop, const G&)
{
  typedef typename boost::graph_traits<G>::face_descriptor face_descriptor;
  return internal::Dynamic_property_map<face_descriptor,T>(prop.t);
}

template<class G, class T, typename Descriptor>
void remove_property(
  internal::Dynamic_property_map<Descriptor, T> pm,
  const G&)
{
  pm.clear();
}



} // namespace internal
} // namespace CGAL

#endif // CGAL_DYNAMIC_PROPERTY_MAP_H
