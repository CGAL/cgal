// Copyright (c) 2017  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_DYNAMIC_PROPERTY_MAP_H
#define CGAL_DYNAMIC_PROPERTY_MAP_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/property_map.h>

#include <memory>

#include <type_traits>
#include <unordered_map>

namespace CGAL {

namespace internal {

template <typename K, typename V>
struct Dynamic_property_map {

  typedef K key_type;
  typedef V value_type;
  typedef const value_type& reference;
  typedef boost::read_write_property_map_tag  category;

  Dynamic_property_map(const V default_value = V())
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


  typedef std::unordered_map<K,V> Map;
  std::shared_ptr<Map> map_;
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
  typedef boost::read_write_property_map_tag category;

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

  std::shared_ptr<PM> map_;
};

template <typename Key, typename Value>
struct Dynamic_with_index
{
  typedef Key key_type;
  typedef Value value_type;
  typedef std::conditional_t<  std::is_same_v<bool, Value>,
                               value_type,
                               value_type&> reference;
  typedef boost::read_write_property_map_tag category;

  Dynamic_with_index()
    : m_values()
  {}

  Dynamic_with_index(std::size_t num_features, Value default_value = Value())
    : m_values( new std::vector<value_type>(num_features, default_value) )
  {}

  friend reference get(const Dynamic_with_index& m, const key_type& k)
  {
    return (*m.m_values)[k.idx()];
  }

  friend void put(const Dynamic_with_index& m, const key_type& k, const value_type& v)
  {
    (*m.m_values)[k.idx()]=v;
  }

  std::shared_ptr<std::vector<value_type> > m_values;
};

} // namespace internal

struct dynamic_property_t {};

template <typename T>
struct dynamic_vertex_property_t : public dynamic_property_t
{
  dynamic_vertex_property_t()
  {}

  using value_type = T;

  template <typename G>
  struct property_map {
    using descriptor = typename boost::graph_traits<G>::vertex_descriptor;
    using vertex_descriptor = descriptor;
  };
};


template <typename T>
struct dynamic_halfedge_property_t : public dynamic_property_t
{
  dynamic_halfedge_property_t()
  {}

  using value_type = T;

  template <typename G>
  struct property_map {
    using descriptor = typename boost::graph_traits<G>::halfedge_descriptor;
    using halfedge_descriptor = descriptor;
  };
};

template <typename T>
struct dynamic_edge_property_t : public dynamic_property_t
{
  dynamic_edge_property_t()
  {}

  using value_type = T;

  template <typename G>
  struct property_map {
    using descriptor = typename boost::graph_traits<G>::edge_descriptor;
    using edge_descriptor = descriptor;
  };
};


template <typename T>
struct dynamic_face_property_t : public dynamic_property_t
{
  dynamic_face_property_t()
  {}

  using value_type = T;

  template <typename G>
  struct property_map
  {
    using descriptor = typename boost::graph_traits<G>::face_descriptor;
    using face_descriptor = descriptor;
  };
};

template <typename T>
constexpr bool is_dynamic_property_tag()
{
  return std::is_base_of_v<CGAL::dynamic_property_t, T>;
}

template <typename Graph, typename Dynamic_property_tag>
struct property_map_of_dynamic_property_map : Dynamic_property_tag::template property_map<Graph>
{
  using descriptor = typename Dynamic_property_tag::template property_map<Graph>::descriptor;
  using type = CGAL::internal::Dynamic_property_map<descriptor, typename Dynamic_property_tag::value_type>;
  using const_type = const type;
};

} // namespace CGAL

namespace boost {

template <typename G, typename T>
struct property_map<G, CGAL::dynamic_vertex_property_t<T> >
  : public CGAL::property_map_of_dynamic_property_map<G, CGAL::dynamic_vertex_property_t<T>>
{
  using vertex_descriptor = typename boost::graph_traits<G>::vertex_descriptor;
};

template <typename G, typename T>
struct property_map<G, CGAL::dynamic_halfedge_property_t<T> >
  : public CGAL::property_map_of_dynamic_property_map<G, CGAL::dynamic_halfedge_property_t<T>>
{
  using halfedge_descriptor = typename boost::graph_traits<G>::halfedge_descriptor;
};


template <typename G, typename T>
struct property_map<G, CGAL::dynamic_edge_property_t<T> >
  : public CGAL::property_map_of_dynamic_property_map<G, CGAL::dynamic_edge_property_t<T>>
{
  using edge_descriptor = typename boost::graph_traits<G>::edge_descriptor;
};

template <typename G, typename T>
struct property_map<G, CGAL::dynamic_face_property_t<T> >
  : public CGAL::property_map_of_dynamic_property_map<G, CGAL::dynamic_face_property_t<T>>
{
  using face_descriptor = typename boost::graph_traits<G>::face_descriptor;
};

} // namespace boost

namespace CGAL {

template <typename Dynamic_property_tag,
          typename G,
          typename = std::enable_if_t<is_dynamic_property_tag<Dynamic_property_tag>()>,
          typename ...Default_value_args>

auto get(const Dynamic_property_tag&, const G&, Default_value_args&&... default_value_args)
{
  using Property_map =
      internal::Dynamic_property_map<typename Dynamic_property_tag::descriptor,
                                     typename Dynamic_property_tag::value_type>;
  return Property_map(std::forward<Default_value_args>(default_value_args)...);
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
