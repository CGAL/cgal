// Copyright (c) 2017  GeometryFactory (France).  All rights reserved.
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
// SPDX-License-Identifier: LGPL-3.0+
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
      (*(const_cast<Dynamic_property_map&>(m).map_))[k] = m.default_value();
      return m.default_value();
    }
    return it->second;
  }


  friend void put(const Dynamic_property_map& m, const key_type& k, const value_type& v)
  {
    (*(m.map_))[k] = v;
  }

  
  const V& default_value() const
  {
    return default_value_;
  }


  typedef boost::unordered_map<K,V> Map;
  boost::shared_ptr<Map> map_;
  V default_value_;
};

  
template <typename M, typename PM>
struct Dynamic_property_map_deleter {
  M& mesh;

  Dynamic_property_map_deleter(const M& mesh)
    : mesh(const_cast<M&>(mesh))
  {}

  void operator()(PM* pm) const
  {
    remove_property(*pm, mesh);
    delete pm;
  }
};


template <typename Mesh, typename PM>
struct Dynamic {
  typedef typename PM::key_type key_type;
  typedef typename PM::value_type value_type;
  typedef typename PM::reference reference;
  typedef typename PM::category category;

  typedef Dynamic_property_map_deleter<Mesh,PM> Deleter;

  Dynamic()
    : map_()
  {}
  
  Dynamic(const Mesh& mesh, PM* pm)
    : map_(pm, Deleter(mesh))
  {}
             
  friend reference get(const Dynamic& m, const key_type& k)
  {
    return get(*(m.map_), k);
  }
    

  friend void put(const Dynamic& m, const key_type& k, const value_type& v)
  {
    put(*(m.map_), k, v);
  }
   
  boost::shared_ptr<PM> map_;
};
  
} // namespace internal

  
template <typename T>
struct dynamic_vertex_property_t
{
  dynamic_vertex_property_t()
  {}
};


template <typename T>
struct dynamic_halfedge_property_t
{
  dynamic_halfedge_property_t()
  {}
};

template <typename T>
struct dynamic_edge_property_t
{
  dynamic_edge_property_t()
  {}
};


template <typename T>
struct dynamic_face_property_t
{
  dynamic_face_property_t()
  {}
};

} // namespace CGAL

namespace boost {


template <typename G, typename T>
struct property_map<G, CGAL::dynamic_vertex_property_t<T> >
{
  typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
  typedef CGAL::internal::Dynamic_property_map<vertex_descriptor,T> type;
  typedef type const_type;
};

template <typename G, typename T>
struct property_map<G, CGAL::dynamic_halfedge_property_t<T> >
{
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<halfedge_descriptor,T> type;
  typedef type const_type;
};


template <typename G, typename T>
struct property_map<G, CGAL::dynamic_edge_property_t<T> >
{
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<edge_descriptor,T> type;
  typedef type const_type;
};

template <typename G, typename T>
struct property_map<G, CGAL::dynamic_face_property_t<T> >
{
  typedef typename boost::graph_traits<G>::face_descriptor face_descriptor;
  typedef CGAL::internal::Dynamic_property_map<face_descriptor,T> type;
  typedef type const_type;
};

} // namespace boost

namespace CGAL {

template <typename T, typename G>
typename boost::property_map<G, dynamic_vertex_property_t<T> >::const_type
get(const CGAL::dynamic_vertex_property_t<T>&, const G&)
{
  typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
  return internal::Dynamic_property_map<vertex_descriptor,T>();
}

template <typename T, typename G>
typename boost::property_map<G, dynamic_halfedge_property_t<T> >::const_type
get(const CGAL::dynamic_halfedge_property_t<T>&, const G&)
{
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  return internal::Dynamic_property_map<halfedge_descriptor,T>();
}

template <typename T, typename G>
typename boost::property_map<G, dynamic_edge_property_t<T> >::const_type
get(const CGAL::dynamic_edge_property_t<T>&, const G&)
{
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  return internal::Dynamic_property_map<edge_descriptor,T>();
}

template <typename T, typename G>
typename boost::property_map<G, dynamic_face_property_t<T> >::const_type
get(const CGAL::dynamic_face_property_t<T>&, const G&)
{
  typedef typename boost::graph_traits<G>::face_descriptor face_descriptor;
  return internal::Dynamic_property_map<face_descriptor,T>();
}

template<typename G, typename Descriptor, typename T>
void remove_property(
  internal::Dynamic_property_map<Descriptor, T> pm,
  const G&)
{
  pm.clear();
}




} // namespace CGAL

#endif // CGAL_DYNAMIC_PROPERTY_MAP_H
