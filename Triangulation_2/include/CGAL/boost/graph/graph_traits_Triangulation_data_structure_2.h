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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_GRAPH_TRAITS_TRIANGULATION_DATA_STRUCTURE_2_H
#define CGAL_GRAPH_TRAITS_TRIANGULATION_DATA_STRUCTURE_2_H

#include <functional>

// include this to avoid a VC15 warning
#include <CGAL/boost/graph/named_function_params.h>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/iterator.h>

// The functions and classes in this file allows the user to
// treat a CGAL Triangulation_data_structure_2 object as a boost graph "as is". No
// wrapper is needed for the Triangulation_data_structure_2 object.

namespace CGAL {

namespace detail {

// A TDS edge is a face handle + an int, and is thus actually a halfedge...
template <class TDS>
struct TDS2_halfedge_descriptor
  : public TDS::Edge
{
  typedef typename TDS::Edge                                   Base;
  typedef typename TDS::Face_handle                            Face_handle;

  TDS2_halfedge_descriptor() {}
  TDS2_halfedge_descriptor(Face_handle fh, int i) : Base(fh, i) { }
  explicit TDS2_halfedge_descriptor(const Base& e) : Base(e) { }
  TDS2_halfedge_descriptor(const TDS2_halfedge_descriptor& h) : Base(h) { }

  const Base& base() const { return static_cast<const Base&>(*this); }

  TDS2_halfedge_descriptor& operator=(const TDS2_halfedge_descriptor& h)
  {
    this->first = h.first;
    this->second = h.second;
    return *this;
  }

  friend std::size_t hash_value(const TDS2_halfedge_descriptor& e) {
    return hash_value(e.first);
  }

  bool operator==(const TDS2_halfedge_descriptor& other) const {
    return (this->first == other.first) && (this->second == other.second);
  }
  bool operator!=(const TDS2_halfedge_descriptor& other) const {
    return (this->first != other.first) || (this->second != other.second);
  }

  bool operator<(const TDS2_halfedge_descriptor& other) const
  {
    if(this->first < other.first) return true;
    if(this->first > other.first) return false;
    return this->second < other.second;
  }
};

// An edge is just a halfedge, but we give it a complete structure to distinguish it from TDS::Edge
template <typename TDS>
struct TDS2_edge_descriptor
{
  typedef typename TDS::Face_handle                              Face_handle;

  TDS2_edge_descriptor() : first(), second(0) { }
  explicit TDS2_edge_descriptor(const typename TDS::Edge& e) : first(e.first), second(e.second) { }
  TDS2_edge_descriptor(Face_handle fd, int i) : first(fd), second(i) { }

  // so that we can still do stuff like tr.is_finite(edge_descriptor) without any hassle
  operator std::pair<Face_handle, int>() { return std::make_pair(first, second); }

  friend std::size_t hash_value(const TDS2_edge_descriptor& h)
  {
    if(h.first == Face_handle())
      return 0;

    return hash_value(h.first < h.first->neighbor(h.second) ? h.first
                                                            : h.first->neighbor(h.second));
  }

  bool operator==(const TDS2_edge_descriptor& other) const
  {
    if((first == other.first) && (second == other.second))
      return true;

    Face_handle fh = first->neighbor(second);
    if(other.first != fh)
      return false;

    int i = fh->index(first);
    return (other.second == i);
  }
  bool operator!=(TDS2_edge_descriptor& other) const { return ! (*this == other); }

  Face_handle first;
  int second;
};

// A halfedge iterator is just an edge iterator that duplicates everything twice,
// to see the edge from either side.
// Could probably be factorized with T2_edge_iterator, but it's clearer this way.
template <typename TDS>
struct TDS2_halfedge_iterator
{
private:
  typedef TDS2_halfedge_iterator<TDS>                           Self;
  typedef typename TDS::Edge_iterator                           Edge_iterator;
  typedef TDS2_halfedge_descriptor<TDS>                         Descriptor;
  typedef typename TDS::Face_handle                             Face_handle;

public:
  typedef Descriptor                                            value_type;
  typedef value_type*                                           pointer;
  typedef value_type&                                           reference;
  typedef std::size_t                                           size_type;
  typedef std::ptrdiff_t                                        difference_type;
  typedef std::bidirectional_iterator_tag                       iterator_category;

  TDS2_halfedge_iterator() { }
  TDS2_halfedge_iterator(const Edge_iterator& feit) : it(feit), on_adjacent_face(false) { }

  Self& operator++()
  {
    // If we are on the first face, move to the opposite face. If we are already on the opposite face,
    // then it's time to move on the next edge
    if(on_adjacent_face) {
      ++it;
      on_adjacent_face = false;
    } else {
      on_adjacent_face = true;
    }

    return *this;
  }

  Self& operator--()
  {
    // Note that while decreasing, we start from the opposite face
    if(on_adjacent_face) {
      on_adjacent_face = false;
    } else {
      --it;
      on_adjacent_face = true;
    }

    return *this;
  }

  Self operator++(int) { Self tmp = *this; operator++(); return tmp; }
  Self operator--(int) { Self tmp = *this; operator--(); return tmp; }

  bool operator==(const Self& other) const { return it == other.it; }
  bool operator!=(const Self& other) const { return !(*this == other); }

  value_type operator*() const
  {
    if(on_adjacent_face) {
      Face_handle neigh_f = it->first->neighbor(it->second);
      return value_type(neigh_f, neigh_f->index(it->first));
     } else {
      return value_type(*it);
    }
  }

private:
  Edge_iterator it;
  bool on_adjacent_face;
};

template <typename TDS>
struct TDS2_edge_iterator
{
private:
  typedef TDS2_edge_iterator<TDS>                               Self;
  typedef typename TDS::Edge_iterator                           Edge_iterator;
  typedef TDS2_edge_descriptor<TDS>                             Descriptor;

public:
  typedef Descriptor                                            value_type;
  typedef value_type*                                           pointer;
  typedef value_type&                                           reference;
  typedef std::size_t                                           size_type;
  typedef std::ptrdiff_t                                        difference_type;
  typedef std::bidirectional_iterator_tag                       iterator_category;

  TDS2_edge_iterator() { }
  TDS2_edge_iterator(const Edge_iterator& feit) : it(feit) { }

  bool operator==(const Self& other) const { return it == other.it; }
  bool operator!=(const Self& other) const { return !(*this == other); }
  Self& operator++() { ++it; return *this; }
  Self& operator--() { --it; return *this; }
  Self operator++(int) { Self tmp = *this; operator++(); return tmp; }
  Self operator--(int) { Self tmp = *this; operator--(); return tmp; }
  value_type operator*() const { return value_type(*it); }

private:
  Edge_iterator it;
};

template <class Circ, class E>
class TDS2_Out_edge_circulator
  : public Circ
{
private:
  mutable E e;

public:
  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  TDS2_Out_edge_circulator() : Circ() {}
  TDS2_Out_edge_circulator(Circ c) : Circ(c) {}

  const E& operator*() const
  {
    E ed = static_cast<const Circ*>(this)->operator*();
    e = E(ed.first->neighbor(ed.second), ed.first->neighbor(ed.second)->index(ed.first));
    return e;
  }
};

template <class Circ, class E>
class TDS2_In_edge_circulator
  : public Circ
{
private:
  mutable E e;

public:
  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  TDS2_In_edge_circulator() : Circ() {}
  TDS2_In_edge_circulator(Circ c) : Circ(c) {}

  const E& operator*() const
  {
    typename Circ::value_type ed = static_cast<const Circ*>(this)->operator*();
    e = E(ed);
    return e;
  }
};

} // namespace detail
} // namespace CGAL

namespace std{

// workaround a bug detected on at least g++ 4.4 where boost::next(Iterator)
// is picked as a candidate for next(h,g)
template <typename TDS>
struct iterator_traits< CGAL::detail::TDS2_halfedge_descriptor<TDS> >
{
  typedef void* iterator_category;
  typedef void* difference_type;
  typedef void* value_type;
  typedef void* reference;
};

} // end of namespace std

namespace boost { 

template <class VB, class FB>
struct graph_traits<CGAL::Triangulation_data_structure_2<VB, FB> >
{
  struct TDS2_graph_traversal_category :
      public virtual bidirectional_graph_tag,
      public virtual adjacency_graph_tag,
      public virtual edge_list_graph_tag,
      public virtual vertex_list_graph_tag { };

  typedef CGAL::Triangulation_data_structure_2<VB,FB> Triangulation_data_structure;

  typedef typename Triangulation_data_structure::Vertex_handle                  vertex_descriptor;
  typedef CGAL::detail::TDS2_halfedge_descriptor<Triangulation_data_structure>  halfedge_descriptor;
  typedef CGAL::detail::TDS2_edge_descriptor<Triangulation_data_structure>      edge_descriptor;
  typedef typename Triangulation_data_structure::Face_handle                    face_descriptor;

  typedef CGAL::Prevent_deref<typename Triangulation_data_structure::Vertex_iterator> vertex_iterator;
  typedef CGAL::detail::TDS2_halfedge_iterator<Triangulation_data_structure>          halfedge_iterator;
  typedef CGAL::detail::TDS2_edge_iterator<Triangulation_data_structure>              edge_iterator;
  typedef CGAL::Prevent_deref<typename Triangulation_data_structure::Face_iterator>   face_iterator;

  typedef CGAL::Counting_iterator<CGAL::detail::TDS2_Out_edge_circulator<typename Triangulation_data_structure::Edge_circulator, edge_descriptor>, edge_descriptor > out_edge_iterator;
  typedef CGAL::Counting_iterator<CGAL::detail::TDS2_In_edge_circulator<typename Triangulation_data_structure::Edge_circulator, edge_descriptor>, edge_descriptor > in_edge_iterator;
  typedef CGAL::Counting_iterator<typename Triangulation_data_structure::Vertex_circulator> Incident_vertices_iterator;
  typedef Incident_vertices_iterator adjacency_iterator;

  typedef undirected_tag directed_category;
  typedef disallow_parallel_edge_tag edge_parallel_category;
  typedef TDS2_graph_traversal_category traversal_category;

  typedef typename Triangulation_data_structure::size_type                          size_type;
  typedef size_type                                                                 vertices_size_type;
  typedef size_type                                                                 edges_size_type;
  typedef size_type                                                                 halfedges_size_type;
  typedef size_type                                                                 faces_size_type;
  typedef size_type                                                                 degree_size_type;

  // nulls
  static vertex_descriptor   null_vertex() { return vertex_descriptor(); }
  static face_descriptor     null_face()   { return face_descriptor(); }
  static halfedge_descriptor     null_halfedge()   { return halfedge_descriptor(); }
};

} // namespace boost

namespace CGAL {

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
next(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
     const Triangulation_data_structure_2<VB,FB>& )
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(e.first, Triangulation_data_structure_2<VB,FB>::ccw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
prev(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
     const Triangulation_data_structure_2<VB,FB>& )
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(e.first, Triangulation_data_structure_2<VB,FB>::cw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
opposite(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
         const Triangulation_data_structure_2<VB,FB>& g)
{
  typedef typename Triangulation_data_structure_2<VB,FB>::Edge Edge;
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(g.mirror_edge(Edge(e.first, e.second)));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor
source(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_descriptor e,
       const Triangulation_data_structure_2<VB,FB>& )
{
  return e.first->vertex(Triangulation_data_structure_2<VB,FB>::ccw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor
target(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_descriptor e,
       const Triangulation_data_structure_2<VB,FB>& )
{
  return e.first->vertex(Triangulation_data_structure_2<VB,FB>::cw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor
source(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
       const Triangulation_data_structure_2<VB,FB>& )
{
  return e.first->vertex(Triangulation_data_structure_2<VB,FB>::ccw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor
target(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
       const Triangulation_data_structure_2<VB,FB>& )
{
  return e.first->vertex(Triangulation_data_structure_2<VB,FB>::cw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::face_descriptor
face(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
     const Triangulation_data_structure_2<VB,FB>&)
{
  return e.first;
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::face_descriptor f,
         const Triangulation_data_structure_2<VB,FB>&)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(f,0);
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor v,
         const Triangulation_data_structure_2<VB,FB>& )
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::face_descriptor face_descriptor;
  face_descriptor fd = v->face();
  int i = fd->index(v);
  return halfedge_descriptor(fd,Triangulation_data_structure_2<VB,FB>::ccw(i));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_descriptor e,
         const Triangulation_data_structure_2<VB,FB>&)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(e.first, e.second);
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_descriptor
edge(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
     const Triangulation_data_structure_2<VB,FB>&)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_descriptor edge_descriptor;
  return edge_descriptor(e.first,e.second);
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_iterator>
vertices(const Triangulation_data_structure_2<VB,FB>& g)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_iterator Iter;
  return make_range(Iter(g.vertices_begin()), Iter(g.vertices_end()));
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_iterator>
edges(const Triangulation_data_structure_2<VB,FB>& g)
{
  typedef typename boost::graph_traits<Triangulation_data_structure_2<VB,FB> >::edge_iterator Iter;
  return make_range(Iter(g.edges_begin()), Iter(g.edges_end()));
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_iterator >
halfedges(const Triangulation_data_structure_2<VB,FB>& g)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_iterator Iter;
  return make_range(Iter(g.edges_begin()), Iter(g.edges_end()));
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::face_iterator >
faces(const Triangulation_data_structure_2<VB,FB>& g)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::face_iterator Iter;
  return make_range(Iter(g.faces_begin()), Iter(g.faces_end()));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type
out_degree(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
           const Triangulation_data_structure_2<VB,FB>& g)
{
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type deg = 0;
  typename Triangulation_data_structure_2<VB,FB>::Edge_circulator c = g.incident_edges(u), done(c);
  if ( c != 0) {
    do {
      ++deg;
    } while (++c != done);
  }
  return deg;
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::out_edge_iterator >
out_edges(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
          const Triangulation_data_structure_2<VB,FB>& g)
{
  typename Triangulation_data_structure_2<VB,FB>::Edge_circulator ec(u,u->face());
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type out_deg = out_degree(u,g);
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::out_edge_iterator Iter;

  return make_range( Iter(ec), Iter(ec,out_deg) );
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::in_edge_iterator >
in_edges(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
         const Triangulation_data_structure_2<VB,FB>& g)
{
  typename Triangulation_data_structure_2<VB,FB>::Edge_circulator ec(u,u->face());
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type out_deg = out_degree(u,g);
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::in_edge_iterator Iter;
  return make_range( Iter(ec), Iter(ec,out_deg) );
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::adjacency_iterator>
adjacent_vertices(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
                  const Triangulation_data_structure_2<VB,FB>& g)
{
  typename Triangulation_data_structure_2<VB,FB>::Vertex_circulator vc = out_edge_iterator(u,u.face());
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type out_deg = out_degree(u,g);
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::adjacency_iterator Iter;
  return make_range( Iter(vc), Iter(vc,out_deg) );
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertices_size_type
num_vertices(const Triangulation_data_structure_2<VB,FB>& g)
{
  return g.number_of_vertices();
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edges_size_type
num_edges(const Triangulation_data_structure_2<VB,FB>& g)
{
  return  g.number_of_vertices() + g.number_of_faces() - 2;
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedges_size_type
num_halfedges(const Triangulation_data_structure_2<VB,FB>& g)
{
  return num_edges(g) * 2;
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::faces_size_type
num_faces(const Triangulation_data_structure_2<VB,FB>& g)
{
  return g.number_of_faces();
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type
in_degree(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
          const Triangulation_data_structure_2<VB,FB>& g)
{
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type deg = 0;
  typename Triangulation_data_structure_2<VB,FB>::Edge_circulator c = g.incident_edges(u), done(c);
  if ( c != 0) {
    do {
      ++deg;
    } while (++c != done);
  }
  return deg;
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type
degree(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
       const Triangulation_data_structure_2<VB,FB>& g)
{
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type deg = 0;
  typename Triangulation_data_structure_2<VB,FB>::Edge_circulator c = g.incident_edges(u), done(c);
  if ( c != 0) {
    do {
      ++deg;
    } while (++c != done);
  }
  return deg;
}

// property maps
template <class VB, class FB>
class TDS2_vertex_id_map
    : public boost::put_get_helper<int, TDS2_vertex_id_map<VB,FB> >
{
public:
  typedef boost::readable_property_map_tag category;
  typedef int value_type;
  typedef int reference;
  typedef typename CGAL::Triangulation_data_structure_2<VB,FB>::Vertex_handle key_type;

  TDS2_vertex_id_map() {}

  long operator[](key_type vh) const { return vh->id(); }
};

template <class VB, class FB>
class TDS2_vertex_point_map
{
public:
  typedef boost::lvalue_property_map_tag category;
  typedef typename VB::Point value_type;
  typedef value_type& reference;
  typedef typename CGAL::Triangulation_data_structure_2<VB,FB>::Vertex_handle key_type;

  friend reference get(TDS2_vertex_point_map<VB,FB>, key_type vh) { return vh->point(); }
  friend void put(TDS2_vertex_point_map<VB,FB>, key_type vh, reference v) { vh->point() = v; }
  reference operator[](key_type vh) const { return vh->point(); }
};

template <class VB, class FB>
class TDS2_edge_id_map
  : public boost::put_get_helper<int, TDS2_edge_id_map<VB,FB> >
{
public:
  typedef boost::readable_property_map_tag category;
  typedef int value_type;
  typedef int reference;
  typedef typename CGAL::Triangulation_data_structure_2<VB,FB>::Edge key_type;

  TDS2_edge_id_map() {}

  long operator[](key_type e) const { return (3 * e.first.id()) + e.second; }
};

template <class VB, class FB>
class TDS2_edge_weight_map
  : public boost::put_get_helper<typename VB::FT, TDS2_edge_weight_map<VB, FB> >
{
private:
  const CGAL::Triangulation_data_structure_2<VB,FB>& tr;

public:
  typedef boost::readable_property_map_tag category;
  typedef typename VB::FT value_type;
  typedef value_type reference;
  typedef typename CGAL::Triangulation_data_structure_2<VB,FB>::Edge key_type;

  TDS2_edge_weight_map(const CGAL::Triangulation_data_structure_2<VB,FB>& tr_) : tr(tr_) { }

  value_type operator[](key_type e) const { return approximate_sqrt(tr.segment(e).squared_length()); }
};

template <class VB, class FB>
inline TDS2_vertex_id_map<VB,FB>
get(boost::vertex_index_t, const Triangulation_data_structure_2<VB,FB>&)
{
  TDS2_vertex_id_map<VB,FB> m;
  return m;
}

template <class VB, class FB>
inline TDS2_vertex_point_map<VB,FB>
get(boost::vertex_point_t, const Triangulation_data_structure_2<VB,FB>&)
{
  TDS2_vertex_point_map<VB,FB> m;
  return m;
}

template <class VB, class FB>
inline TDS2_edge_id_map<VB,FB>
get(boost::edge_index_t, const Triangulation_data_structure_2<VB,FB>&)
{
  TDS2_edge_id_map<VB,FB> m;
  return m;
}

template <class VB, class FB>
inline TDS2_edge_weight_map<VB,FB>
get(boost::edge_weight_t, const Triangulation_data_structure_2<VB,FB>& g)
{
  TDS2_edge_weight_map<VB,FB> m(g);
  return m;
}

template <class Tag>
struct TDS2_property_map { };

template <>
struct TDS2_property_map<boost::vertex_index_t>
{
  template <class VB, class FB>
  struct bind_ {
    typedef TDS2_vertex_id_map<VB,FB> type;
    typedef TDS2_vertex_id_map<VB,FB> const_type;
  };
};



template <>
struct TDS2_property_map<boost::vertex_point_t>
{
  template <class VB, class FB>
  struct bind_ {
    typedef TDS2_vertex_point_map<VB,FB> type;
    typedef TDS2_vertex_point_map<VB,FB> const_type;
  };
};


template <>
struct TDS2_property_map<boost::edge_index_t>
{
  template <class VB, class FB>
  struct bind_ {
    typedef TDS2_edge_id_map<VB,FB> type;
    typedef TDS2_edge_id_map<VB,FB> const_type;
  };
};


template <>
struct TDS2_property_map<boost::edge_weight_t>
{
  template <class VB, class FB>
  struct bind_ {
    typedef TDS2_edge_weight_map<VB,FB> type;
    typedef TDS2_edge_weight_map<VB,FB> const_type;
  };
};

} // namespace CGAL

namespace boost {
// g++ 'enumeral_type' in template unification not implemented workaround
template <class VB, class FB, class Tag>
struct property_map<CGAL::Triangulation_data_structure_2<VB,FB>, Tag>
{
  typedef typename
  CGAL::TDS2_property_map<Tag>::template bind_<VB,FB> map_gen;
  typedef typename map_gen::type type;
  typedef typename map_gen::const_type const_type;
};

// see struct property_map in Polyehdron for an explanation
template <class VB, class FB, class Tag>
struct property_map<const CGAL::Triangulation_data_structure_2<VB,FB>, Tag>
{
  typedef typename
  CGAL::TDS2_property_map<Tag>::template bind_<VB,FB> map_gen;
  typedef typename map_gen::type type;
  typedef typename map_gen::const_type const_type;
};

} // namespace boost


namespace CGAL {

template <class VB, class FB, class PropertyTag, class Key>
inline
typename boost::property_traits<
typename boost::property_map<Triangulation_data_structure_2<VB,FB>,PropertyTag>::const_type>::value_type
get(PropertyTag p, const Triangulation_data_structure_2<VB,FB>& g, const Key& key)
{
  return get(get(p, g), key);
}

template <class VB, class FB, class PropertyTag, class Key,class Value>
inline void
put(PropertyTag p, Triangulation_data_structure_2<VB,FB>& g,
    const Key& key, const Value& value)
{
  typedef typename boost::property_map<Triangulation_data_structure_2<VB,FB>, PropertyTag>::type Map;
  Map pmap = get(p, g);
  put(pmap, key, value);
}

} // namespace CGAL 

namespace boost {

// What are those needed for ???
template <typename VB, typename FB>
struct edge_property_type<CGAL::Triangulation_data_structure_2<VB,FB> > {
  typedef void type;
};

template <typename VB, typename FB>
struct vertex_property_type<CGAL::Triangulation_data_structure_2<VB,FB> > {
  typedef void type;
};
} // namespace boost

namespace std {

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4099) // For VC10 it is class hash 
#endif

#ifndef CGAL_CFG_NO_STD_HASH

template <typename TDS>
struct hash<CGAL::detail::TDS2_edge_descriptor<TDS> >
{
  std::size_t operator()(const CGAL::detail::TDS2_edge_descriptor<TDS>& e) const { return hash_value(e); }
};

template <typename TDS>
struct hash<CGAL::detail::TDS2_halfedge_descriptor<TDS> >
{
  std::size_t operator()(const CGAL::detail::TDS2_halfedge_descriptor<TDS>& e) const { return hash_value(e); }
};

#endif // CGAL_CFG_NO_STD_HASH

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // namespace std

#endif // CGAL_GRAPH_TRAITS_TRIANGULATION_DATA_STRUCTURE_2_H
