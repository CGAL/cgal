// Copyright (c) 2007, 2019  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS
  #error CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS is not defined
#endif

#ifndef CGAL_2D_TRIANGULATION
  #error CGAL_2D_TRIANGULATION is not defined
#endif

#ifndef CGAL_2D_TRIANGULATION_TEMPLATES
  #error CGAL_2D_TRIANGULATION_TEMPLATES is not defined
#endif

#include <functional>

// include this to avoid a VC15 warning
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_data_structure_2.h>
#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/iterator.h>

// Small include guard is required for that class
#ifndef CGAL_GRAPH_TRAITS_2D_TRIANGULATION
#define CGAL_GRAPH_TRAITS_2D_TRIANGULATION

namespace CGAL {

namespace detail {

template <typename Tr>
struct T2_halfedge_descriptor
{
  typedef typename Tr::Face_handle face_descriptor;
  face_descriptor first;
  int second;
  operator std::pair<face_descriptor, int>() { return std::make_pair(first, second); }

  T2_halfedge_descriptor() : first(), second(0) { }
  T2_halfedge_descriptor(const typename Tr::Edge& e) : first(e.first), second(e.second) { }
  T2_halfedge_descriptor(face_descriptor fd, int i) : first(fd), second(i) { }

  friend std::size_t hash_value(const T2_halfedge_descriptor& h)
  {
    return hash_value(h.first);
  }

  bool operator==(const T2_halfedge_descriptor& other) const
  {
    return (first == other.first) && (second == other.second);
  }

  bool operator!=(const T2_halfedge_descriptor& other) const
  {
    return (first != other.first) || (second != other.second);
  }

  bool operator<(const T2_halfedge_descriptor& other) const
  {
    if(first < other.first) return true;
    if(first > other.first) return false;
    return second  < other.second;
  }
};

template <typename Tr>
class T2_vertex_id_map
  : public boost::put_get_helper<int, T2_vertex_id_map<Tr> >
{
public:
  typedef boost::readable_property_map_tag                         category;
  typedef int                                                      value_type;
  typedef int                                                      reference;
  typedef typename Tr::Vertex_handle                               key_type;

  T2_vertex_id_map() { }

  long operator[](key_type vh) const { return vh->id(); }
};

template <typename Tr>
class T2_vertex_point_map
{
public:
  typedef boost::lvalue_property_map_tag                           category;
  typedef typename Tr::Point                                       value_type;
  typedef value_type&                                              reference;
  typedef typename Tr::Vertex_handle                               key_type;

  friend reference get(T2_vertex_point_map<Tr>, key_type vh)
  {
    return vh->point();
  }
  friend void put(T2_vertex_point_map<Tr>, key_type vh, reference v)
  {
    vh->point() = v;
  }
  reference operator[](key_type vh) const { return vh->point(); }
};

template <typename Tr>
class T2_edge_id_map
  : public boost::put_get_helper<int, T2_edge_id_map<Tr> >
{
public:
  typedef boost::readable_property_map_tag                         category;
  typedef int                                                      value_type;
  typedef int                                                      reference;
  typedef typename Tr::Edge                                        key_type;

  T2_edge_id_map() { }

  long operator[](key_type e) const { return (3 * e.first.id()) + e.second; }
};

template <typename Tr>
class T2_edge_weight_map
  : public boost::put_get_helper<typename Tr::Geom_traits::FT,
                                 T2_edge_weight_map<Tr> >
{
private:
  const Tr& tr;

public:
  typedef boost::readable_property_map_tag                        category;
  typedef typename Tr::Geom_traits::FT                            value_type;
  typedef value_type                                              reference;
  typedef typename Tr::Edge                                       key_type;

  T2_edge_weight_map(const Tr& tr_) : tr(tr_) { }

  value_type operator[](key_type e) const { return approximate_sqrt(tr.segment(e).squared_length()); }
};

template <typename Tr, typename Tag>
struct T2_property_map { };

template <typename Tr>
struct T2_property_map<Tr, boost::vertex_index_t>
{
  typedef detail::T2_vertex_id_map<Tr> type;
  typedef detail::T2_vertex_id_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::vertex_point_t>
{
  typedef detail::T2_vertex_point_map<Tr> type;
  typedef detail::T2_vertex_point_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::edge_index_t>
{
  typedef detail::T2_edge_id_map<Tr> type;
  typedef detail::T2_edge_id_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::edge_weight_t>
{
  typedef detail::T2_edge_weight_map<Tr> type;
  typedef detail::T2_edge_weight_map<Tr> const_type;
};

} // namespace detail
} // namespace CGAL

namespace std {

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4099) // For VC10 it is class hash
#endif

#ifndef CGAL_CFG_NO_STD_HASH

template < class Tr>
struct hash<CGAL::detail::T2_halfedge_descriptor<Tr> >
{
  std::size_t operator()(const CGAL::detail::T2_halfedge_descriptor<Tr>& e) const
  {
    return hash_value(e);
  }
};

#endif // CGAL_CFG_NO_STD_HASH

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // namespace std

#endif // CGAL_GRAPH_TRAITS_2D_TRIANGULATION

namespace boost {

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct graph_traits< CGAL_2D_TRIANGULATION >
{
  typedef CGAL_2D_TRIANGULATION                                               Triangulation;

  struct T2_graph_traversal_category :
      public virtual bidirectional_graph_tag,
      public virtual adjacency_graph_tag,
      public virtual edge_list_graph_tag,
      public virtual vertex_list_graph_tag { };

  typedef typename Triangulation::Vertex_handle                               vertex_descriptor;
  typedef CGAL::detail::T2_halfedge_descriptor<Triangulation>                 halfedge_descriptor;
  typedef CGAL::detail::Edge<Triangulation, typename Triangulation::Edge>     edge_descriptor;
  typedef typename Triangulation::Face_handle                                 face_descriptor;

  typedef CGAL::Prevent_deref<typename Triangulation::All_vertices_iterator>  vertex_iterator;
  typedef typename Triangulation::All_halfedges_iterator                      halfedge_iterator;
  typedef typename Triangulation::All_edges_iterator                          edge_iterator;
  typedef CGAL::Prevent_deref<typename Triangulation::All_faces_iterator>     face_iterator;

  typedef CGAL::Counting_iterator<
  CGAL::detail::Out_edge_circulator<
  typename Triangulation::Edge_circulator,
  edge_descriptor>,
  edge_descriptor>                                                            out_edge_iterator;
  typedef CGAL::Counting_iterator<
  CGAL::detail::In_edge_circulator<
  typename Triangulation::Edge_circulator,
  edge_descriptor>,
  edge_descriptor>                                                            in_edge_iterator;

  typedef CGAL::Counting_iterator<typename Triangulation::Vertex_circulator>  Incident_vertices_iterator;
  typedef Incident_vertices_iterator                                          adjacency_iterator;

  typedef undirected_tag                                                      directed_category;
  typedef disallow_parallel_edge_tag                                          edge_parallel_category;
  typedef T2_graph_traversal_category                                         traversal_category;

  typedef typename Triangulation::size_type                                   size_type;
  typedef size_type                                                           vertices_size_type;
  typedef size_type                                                           halfedges_size_type;
  typedef size_type                                                           edges_size_type;
  typedef size_type                                                           faces_size_type;
  typedef size_type                                                           degree_size_type;

  // nulls
  static vertex_descriptor null_vertex() { return vertex_descriptor(); }
  static halfedge_descriptor null_halfedge() { return halfedge_descriptor(); }
  static face_descriptor null_face() { return face_descriptor(); }
};

} // namespace boost

namespace CGAL {

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
next(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor e,
     const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(e.first, g.ccw(e.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
prev(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor e,
     const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(e.first, g.cw(e.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
opposite(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor e,
         const CGAL_2D_TRIANGULATION& g)
{
  return g.mirror_edge(e);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor
source(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor e,
       const CGAL_2D_TRIANGULATION& g)
{
  return e.first->vertex(g.ccw(e.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor
target(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor e,
       const CGAL_2D_TRIANGULATION& g)
{
  return e.first->vertex(g.cw(e.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor
source(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor e,
       const CGAL_2D_TRIANGULATION& g)
{
  return e.first->vertex(g.ccw(e.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor
target(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor e,
       const CGAL_2D_TRIANGULATION& g)
{
  return e.first->vertex(g.cw(e.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor
face(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor e,
     const CGAL_2D_TRIANGULATION&)
{
  return e.first;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
halfedge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor f,
         const CGAL_2D_TRIANGULATION&)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(f,0);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
halfedge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
         const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor face_descriptor;
  face_descriptor fd = v->face();
  int i = fd->index(v);
  return halfedge_descriptor(fd, g.ccw(i));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
halfedge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor e,
         const CGAL_2D_TRIANGULATION&)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(e.first, e.second);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor
edge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor e,
     const CGAL_2D_TRIANGULATION&)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor edge_descriptor;
  return edge_descriptor(e.first, e.second);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
std::pair<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor,
bool>
edge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor u,
     typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
     const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor edge_descriptor;

  typename CGAL_2D_TRIANGULATION::Edge_circulator c = g.incident_edges(u), done(c);
  if(c != 0) {
    do {
      // find the index of the other vertex of *c
      int indv = 3 - c->first->index(u) - c->second;
      if(c->first->vertex(indv) == v)
        return std::make_pair(edge_descriptor(c->first, c->second), true);
    } while (++c != done);
  }

  return std::make_pair(edge_descriptor(), false);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
std::pair<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor,
bool>
halfedge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor u,
         typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
         const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor face_descriptor;

  std::pair<edge_descriptor, bool> eb = edge(u, v, g);

  if(!eb.second)
    return std::make_pair(halfedge_descriptor(), false);

  const edge_descriptor& e = eb.first;

  if(e.first->vertex(g.ccw(e.first->index(u))) == v)
  {
    return std::make_pair(halfedge_descriptor(e.first, e.second), true);
  }
  else
  {
    face_descriptor nf = e.first->neighbor(e.second);
    int idx = nf->index(e.first);
    return std::make_pair(halfedge_descriptor(nf, idx), true);
  }
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_iterator>
vertices(const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_iterator
      Iter;
  return make_range( Iter(g.all_vertices_begin()), Iter(g.all_vertices_end()) );
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_iterator>
edges(const CGAL_2D_TRIANGULATION& g)
{
  return make_range(g.all_edges_begin(), g.all_edges_end());
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_iterator >
halfedges(const CGAL_2D_TRIANGULATION& g)
{
  return make_range(g.all_halfedges_begin(), g.all_halfedges_end());
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_iterator >
faces(const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_iterator
      Iter;
  return make_range( Iter(g.all_faces_begin()), Iter(g.all_faces_end()) );
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type
out_degree(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor u,
           const CGAL_2D_TRIANGULATION& g)
{
  typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type deg = 0;
  typename CGAL_2D_TRIANGULATION::Edge_circulator c = g.incident_edges(u), done(c);
  if(c != 0)
  {
    do {
      ++deg;
    } while (++c != done);
  }
  return deg;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::out_edge_iterator >
out_edges(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor u,
          const CGAL_2D_TRIANGULATION& g)
{
  typename CGAL_2D_TRIANGULATION::Edge_circulator ec(u, u->face());
  typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type out_deg = out_degree(u, g);
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >
      ::out_edge_iterator Iter;

  return make_range(Iter(ec), Iter(ec, out_deg));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::in_edge_iterator >
in_edges(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor u,
         const CGAL_2D_TRIANGULATION& g)
{
  typename CGAL_2D_TRIANGULATION::Edge_circulator ec(u, u->face());
  typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type out_deg = out_degree(u, g);
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >
      ::in_edge_iterator Iter;
  return make_range(Iter(ec), Iter(ec,out_deg));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::adjacency_iterator>
adjacent_vertices(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor u,
                  const CGAL_2D_TRIANGULATION& g)
{
  typename CGAL_2D_TRIANGULATION::Vertex_circulator vc = out_edge_iterator(u, u.face());
  typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type out_deg = out_degree(u, g);
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::adjacency_iterator Iter;
  return make_range( Iter(vc), Iter(vc,out_deg) );
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertices_size_type
num_vertices(const CGAL_2D_TRIANGULATION& g)
{
  return g.tds().number_of_vertices();
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edges_size_type
num_edges(const CGAL_2D_TRIANGULATION& g)
{
  return g.tds().number_of_vertices() + g.tds().number_of_faces() - 2;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedges_size_type
num_halfedges(const CGAL_2D_TRIANGULATION& g)
{
  return num_edges(g) * 2;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::faces_size_type
num_faces(const CGAL_2D_TRIANGULATION& g)
{
  return g.tds().number_of_faces();
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type
in_degree(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor u,
          const CGAL_2D_TRIANGULATION& g)
{
  typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type deg = 0;
  typename CGAL_2D_TRIANGULATION::Edge_circulator c = g.incident_edges(u), done(c);
  if(c != 0)
  {
    do {
      ++deg;
    } while (++c != done);
  }
  return deg;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type
degree(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor u,
       const CGAL_2D_TRIANGULATION& g)
{
  typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type deg = 0;
  typename CGAL_2D_TRIANGULATION::Edge_circulator c = g.incident_edges(u), done(c);
  if(c != 0)
  {
    do {
      ++deg;
    } while (++c != done);
  }
  return deg;
}

// property maps
template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline detail::T2_vertex_id_map< CGAL_2D_TRIANGULATION >
get(boost::vertex_index_t, const CGAL_2D_TRIANGULATION&)
{
  detail::T2_vertex_id_map< CGAL_2D_TRIANGULATION > m;
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline detail::T2_vertex_point_map< CGAL_2D_TRIANGULATION >
get(boost::vertex_point_t, const CGAL_2D_TRIANGULATION&)
{
  detail::T2_vertex_point_map< CGAL_2D_TRIANGULATION > m;
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline detail::T2_edge_id_map< CGAL_2D_TRIANGULATION >
get(boost::edge_index_t, const CGAL_2D_TRIANGULATION&)
{
  detail::T2_edge_id_map< CGAL_2D_TRIANGULATION > m;
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline detail::T2_edge_weight_map< CGAL_2D_TRIANGULATION >
get(boost::edge_weight_t, const CGAL_2D_TRIANGULATION& g)
{
  detail::T2_edge_weight_map< CGAL_2D_TRIANGULATION > m(g);
  return m;
}

} // namespace CGAL

namespace boost {

// g++ 'enumeral_type' in template unification not implemented workaround
template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS, class Tag>
struct property_map<CGAL_2D_TRIANGULATION, Tag>
{
  typedef typename CGAL::detail::T2_property_map<CGAL_2D_TRIANGULATION, Tag>    map_gen;
  typedef typename map_gen::type                                                type;
  typedef typename map_gen::const_type                                          const_type;
};

// see struct property_map in Polyehdron for an explanation
template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS, class Tag>
struct property_map<const CGAL_2D_TRIANGULATION, Tag>
{
  typedef typename CGAL::detail::T2_property_map<CGAL_2D_TRIANGULATION, Tag>    map_gen;
  typedef typename map_gen::type                                                type;
  typedef typename map_gen::const_type                                          const_type;
};

} // namespace boost

namespace CGAL {

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS, class PropertyTag, class Key>
inline
typename boost::property_traits<
typename boost::property_map< CGAL_2D_TRIANGULATION, PropertyTag>::const_type>::value_type
get(PropertyTag p, const CGAL_2D_TRIANGULATION& g, const Key& key)
{
  return get(get(p, g), key);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS, class PropertyTag, class Key, class Value>
inline void
put(PropertyTag p, CGAL_2D_TRIANGULATION& g,
    const Key& key, const Value& value)
{
  typedef typename boost::property_map<CGAL_2D_TRIANGULATION, PropertyTag>::type Map;
  Map pmap = get(p, g);
  put(pmap, key, value);
}

} // namespace CGAL

namespace boost {

// What are those needed for ???
template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct edge_property_type<CGAL_2D_TRIANGULATION >
{
  typedef void type;
};

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct vertex_property_type<CGAL_2D_TRIANGULATION >
{
  typedef void type;
};

} // namespace boost

#undef CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS
#undef CGAL_2D_TRIANGULATION
#undef CGAL_2D_TRIANGULATION_TEMPLATES
