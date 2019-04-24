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
// Author(s)     : Mael Rouxel-Labbé,
//                 Andreas Fabri,
//                 Fernando Cacciola

#ifndef CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS
  #error CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS is not defined
#endif

#ifndef CGAL_2D_TRIANGULATION
  #error CGAL_2D_TRIANGULATION is not defined
#endif

#ifndef CGAL_2D_TRIANGULATION_TEMPLATES
  #error CGAL_2D_TRIANGULATION_TEMPLATES is not defined
#endif

#include <CGAL/Iterator_range.h>
#include <CGAL/iterator.h>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>

#include <utility>

// Small include guard for helper classes
#ifndef CGAL_GRAPH_TRAITS_2D_TRIANGULATION
#define CGAL_GRAPH_TRAITS_2D_TRIANGULATION

namespace CGAL {
namespace detail {

template <typename Tr>
struct T2_halfedge_descriptor
{
  typedef typename Tr::Face_handle                              face_descriptor;

  T2_halfedge_descriptor() : first(), second(0) { }
  explicit T2_halfedge_descriptor(const typename Tr::Edge& e) : first(e.first), second(e.second) { }
  T2_halfedge_descriptor(face_descriptor fd, int i) : first(fd), second(i) { }

  operator std::pair<face_descriptor, int>() { return std::make_pair(first, second); }

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

  face_descriptor first;
  int second;
};

template <class Tr>
struct T2_edge_descriptor
  : public Tr::Edge
{
  typedef typename Tr::Edge                                   Base;
  typedef typename Tr::Face_handle                            Face_handle;

  T2_edge_descriptor() {}
  T2_edge_descriptor(Face_handle fh, int i) : Base(fh, i) { }
  explicit T2_edge_descriptor(const Base& e) : Base(e) { }
  T2_edge_descriptor(const T2_edge_descriptor& e) : Base(e) { }

  T2_edge_descriptor& operator=(const T2_edge_descriptor& e)
  {
    this->first = e.first;
    this->second = e.second;
    return *this;
  }

  friend std::size_t hash_value(const T2_edge_descriptor& e)
  {
    if (e.first==Face_handle()) return 0;
    return hash_value(e.first<e.first->neighbor(e.second)?
                      e.first:e.first->neighbor(e.second));
  }

  bool operator==(const T2_edge_descriptor& other) const
  {
    if((this->first == other.first)&&(this->second == other.second)) return true;
    Face_handle fh = this->first->neighbor(this->second);
    if(other.first != fh) return false;
    int i = fh->index(this->first);
    return (other.second == i);
  }

  bool operator!=(T2_edge_descriptor& other) const
  {
    return ! (*this == other);
  }
};

template <typename Tr, typename Triangulation_iterator, typename Descriptor>
struct T2_iterator
{
private:
  typedef T2_iterator<Tr, Triangulation_iterator, Descriptor>   Self;
  typedef typename Tr::Triangulation_data_structure             Tds;

public:
  typedef Descriptor                                            value_type;
  typedef value_type*                                           pointer;
  typedef value_type&                                           reference;
  typedef std::size_t                                           size_type;
  typedef std::ptrdiff_t                                        difference_type;
  typedef std::bidirectional_iterator_tag                       iterator_category;

  T2_iterator() { }
  T2_iterator(const Triangulation_iterator& it) : it(it) { }

  bool operator==(const Self& other) const { return it == other.it; }
  bool operator!=(const Self& other) const { return !(*this == other);}
  Self& operator++() { ++it; return *this; }
  Self& operator--() { --it; return *this; }
  Self operator++(int) { Self tmp = *this; operator++(); return tmp; }
  Self operator--(int) { Self tmp = *this; operator--(); return tmp; }
  value_type operator*() const { return value_type(*it); }

  Triangulation_iterator it;
};

template <class Circ, class E>
struct Out_edge_circulator
  : public Circ
{
private:
  mutable E e;

public:
  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  Out_edge_circulator() : Circ() { }
  Out_edge_circulator(Circ c) : Circ(c) { }

  const E& operator*() const
  {
    E ed(static_cast<const Circ*>(this)->operator*());
    e = E(ed.first->neighbor(ed.second), ed.first->neighbor(ed.second)->index(ed.first));
    return e;
  }
};

template <class Circ, class E>
struct In_edge_circulator
  : public Circ
{
private:
  mutable E e;

public:
  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  In_edge_circulator() : Circ() { }
  In_edge_circulator(Circ c) : Circ(c) { }

  const E& operator*() const
  {
    e = E(static_cast<const Circ*>(this)->operator*());
    return e;
  }
};

template <typename Tr>
struct Dereference_to_vertex_handle_enforcer
  : public boost::iterator_adaptor<
      Dereference_to_vertex_handle_enforcer<Tr>,
      typename Tr::All_vertices_iterator /*base*/,
      typename Tr::Vertex_handle /*value*/,
      boost::use_default,
      typename Tr::Vertex_handle /*reference*/
    >
{
public:
  typedef typename Tr::Vertex_handle                                                   value_type;

private:
  typedef Dereference_to_vertex_handle_enforcer<Tr>                                    Self;
  typedef typename Tr::All_vertices_iterator                                           I;
  typedef boost::iterator_adaptor<Self, I, value_type, boost::use_default, value_type> Base;

public:
  Dereference_to_vertex_handle_enforcer() { }
  explicit Dereference_to_vertex_handle_enforcer(const I& i) : Base(i) { }

private:
  friend class boost::iterator_core_access;
  value_type dereference() const { return value_type(this->base()); }
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
      public virtual vertex_list_graph_tag
  { };

  typedef typename Triangulation::Vertex_handle                               vertex_descriptor;
  typedef CGAL::detail::T2_halfedge_descriptor<Triangulation>                 halfedge_descriptor;
  typedef CGAL::detail::T2_edge_descriptor<Triangulation>                     edge_descriptor;
  typedef typename Triangulation::Face_handle                                 face_descriptor;

  // Regular_triangulation_2 unfortunately overrides the type 'All_vertices_iterator' due to hidden
  // points, so we can't simply use 'CGAL::Prevent_deref' since 'typeid(*vertex_iterator)' would
  // not be the vertex_descripor type, which creates all types (!) of problems.
  typedef CGAL::detail::Dereference_to_vertex_handle_enforcer<Triangulation>  vertex_iterator;

  // Since 'All_halfedges_...' and 'All_edges_...' have the same value type, we have to wrap them
  // around to enforce the value_type being the correct descriptor
  typedef CGAL::detail::T2_iterator<Triangulation,
                                    typename Triangulation::All_halfedges_iterator,
                                    halfedge_descriptor>                      halfedge_iterator;
  typedef CGAL::detail::T2_iterator<Triangulation,
                                    typename Triangulation::All_edges_iterator,
                                    edge_descriptor>                          edge_iterator;
  typedef CGAL::Prevent_deref<typename Triangulation::All_faces_iterator>     face_iterator;

  typedef CGAL::Counting_iterator<
            CGAL::detail::In_edge_circulator<
              typename Triangulation::Edge_circulator,
              edge_descriptor>,
            edge_descriptor>                                                  in_edge_iterator;
  typedef CGAL::Counting_iterator<
            CGAL::detail::Out_edge_circulator<
              typename Triangulation::Edge_circulator,
              edge_descriptor>,
            edge_descriptor>                                                  out_edge_iterator;

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

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct graph_traits<const CGAL_2D_TRIANGULATION >
  : public graph_traits< CGAL_2D_TRIANGULATION >
{ };

} // namespace boost

namespace CGAL {

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor
source(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor e,
       const CGAL_2D_TRIANGULATION& g)
{
  return e.first->vertex(g.ccw(e.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor
source(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
       const CGAL_2D_TRIANGULATION& g)
{
  return h.first->vertex(g.ccw(h.second));
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
target(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
       const CGAL_2D_TRIANGULATION& g)
{
  return h.first->vertex(g.cw(h.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
next(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
     const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(h.first, g.ccw(h.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
prev(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
     const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(h.first, g.cw(h.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
opposite(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
         const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor     edge_descriptor;

  return halfedge_descriptor(g.mirror_edge(edge_descriptor(h.first, h.second)));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor
face(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
     const CGAL_2D_TRIANGULATION&)
{
  return h.first;
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
edge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
     const CGAL_2D_TRIANGULATION&)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor edge_descriptor;
  return edge_descriptor(h.first, h.second);
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
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_iterator Iter;
  return make_range( Iter(g.all_vertices_begin()), Iter(g.all_vertices_end()) );
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_iterator >
halfedges(const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_iterator Iter;
  return make_range(Iter(g.all_halfedges_begin()), Iter(g.all_halfedges_end()));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_iterator>
edges(const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_iterator Iter;
  return make_range(Iter(g.all_edges_begin()), Iter(g.all_edges_end()));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_iterator >
faces(const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_iterator Iter;
  return make_range(Iter(g.all_faces_begin()), Iter(g.all_faces_end()));
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
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::in_edge_iterator >
in_edges(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor u,
         const CGAL_2D_TRIANGULATION& g)
{
  typename CGAL_2D_TRIANGULATION::Edge_circulator ec(u, u->face());
  typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type out_deg = out_degree(u, g);
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::in_edge_iterator     Iter;
  return make_range(Iter(ec), Iter(ec,out_deg));
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
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::out_edge_iterator Iter;

  return make_range(Iter(ec), Iter(ec, out_deg));
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

} // namespace CGAL

#undef CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS
#undef CGAL_2D_TRIANGULATION
#undef CGAL_2D_TRIANGULATION_TEMPLATES
