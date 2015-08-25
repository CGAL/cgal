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

#ifndef CGAL_GRAPH_TRAITS_TRIANGULATION_2_H
#define CGAL_GRAPH_TRAITS_TRIANGULATION_2_H

#include <functional>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Iterator_range.h>

// The functions and classes in this file allows the user to
// treat a CGAL Triangulation_2 object as a boost graph "as is". No
// wrapper is needed for the Triangulation_2 object.


namespace CGAL {

  namespace detail {

template < class T, class EdgeBase >
class Edge : public EdgeBase {

public:
  typedef typename T::Face_handle Face_handle ;
  
  Edge()
  {}

  Edge(Face_handle  fh, int i)
    : EdgeBase(fh,i)
  {}
  
  Edge(const EdgeBase& e)
    : EdgeBase(e)
  {}

  Edge(const Edge& e)
    : EdgeBase(e)
  {}

  Edge&
  operator=(const Edge& e)
  {
    this->first = e.first;
    this->second = e.second;
    return *this;
  }

  friend std::size_t hash_value(const Edge& e)
  {
    return hash_value(e.first);
  }

  bool operator==(const Edge& other) const
  {
    if((this->first == other.first)&&(this->second == other.second)) return true;
    Face_handle fh = this->first->neighbor(this->second);
    if(other.first != fh) return false;
    int i = fh->index(this->first);
    return (other.second == i);
  }

  bool operator!=(Edge& other) const
  {
    return ! (*this == other);
  }
};

template <class Circ, class E>
class Out_edge_circulator : public Circ
{
private:
  mutable E e;

public:

  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  Out_edge_circulator()
    : Circ()
  {}

  Out_edge_circulator(Circ c)
    : Circ(c)
  {}

  const E& operator*() const
  {
    E ed = static_cast<const Circ*>(this)->operator*();
    e = E(ed.first->neighbor(ed.second), ed.first->neighbor(ed.second)->index(ed.first));
    return e;
  }
};
 
template <class Circ, class E>
class In_edge_circulator : public Circ
{
private:
  mutable E e;

public:

  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  In_edge_circulator()
    : Circ()
  {}

  In_edge_circulator(Circ c)
    : Circ(c)
  {}

  const E& operator*() const
{
    typename Circ::value_type ed = static_cast<const Circ*>(this)->operator*();
    e = E(ed);
    return e;
  }
};
  

  //  The vertex iterator of the bgl must evaluate to a vertex handle, not to a vertex
template < class T>
class boost_all_vertices_iterator {
protected:
 typename T::All_vertices_iterator nt;
public:
  typedef typename T::All_vertices_iterator  Iterator;
  typedef boost_all_vertices_iterator<T> Self;

  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;
  typedef typename T::Vertex_handle  value_type;
  typedef typename std::iterator_traits<Iterator>::difference_type           difference_type;
  typedef value_type      reference;
  typedef value_type      pointer;

  // CREATION
  // --------

  boost_all_vertices_iterator()
  {}

  boost_all_vertices_iterator( Iterator j) : nt(j) {}

  // OPERATIONS Forward Category
  // ---------------------------


  bool operator==( const Self& i) const { return ( nt == i.nt); }
  bool operator!=( const Self& i) const { return !(nt == i.nt );   }
  value_type  operator*() const  { return nt; }
  value_type    operator->()  { return nt; }

  Self& operator++() {
    ++nt;
    return *this;
  }

  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  Self& operator--() {
    --nt;
    return *this;
  }

  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }
};

template < class T>
class boost_all_faces_iterator {
protected:
 typename T::All_faces_iterator nt;
public:
  typedef typename T::All_faces_iterator  Iterator;
  typedef boost_all_faces_iterator<T> Self;

  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;
  typedef typename T::Face_handle  value_type;
  typedef typename std::iterator_traits<Iterator>::difference_type           difference_type;
  typedef value_type      reference;
  typedef value_type      pointer;

  // CREATION
  // --------

  boost_all_faces_iterator()
  {}

  boost_all_faces_iterator( Iterator j) : nt(j) {}

  // OPERATIONS Forward Category
  // ---------------------------


  bool operator==( const Self& i) const { return ( nt == i.nt); }
  bool operator!=( const Self& i) const { return !(nt == i.nt );   }
  value_type  operator*() const  { return nt; }
  value_type    operator->()  { return nt; }

  Self& operator++() {
    ++nt;
    return *this;
  }

  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  Self& operator--() {
    --nt;
    return *this;
  }

  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }
};

    
    template <typename Tr>
    struct T2_halfedge_descriptor
    {
      typedef typename Tr::Face_handle face_descriptor;
      face_descriptor first;
      int second;
      operator std::pair<face_descriptor, int>() { return std::make_pair(first,second); }
      
      T2_halfedge_descriptor()
      {}
      
      T2_halfedge_descriptor(const typename Tr::Edge& e)
        : first(e.first), second(e.second)
      {}
      
      T2_halfedge_descriptor(face_descriptor fd, int i)
        : first(fd), second(i)
      {}
      
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


  } // namespace detail
} // namespace CGAL

namespace boost { 

  template <class GT, class TDS>
  struct graph_traits< CGAL::Triangulation_2<GT,TDS> > {

    struct T2_graph_traversal_category : 
      public virtual bidirectional_graph_tag,
      public virtual adjacency_graph_tag,
      public virtual edge_list_graph_tag,
      public virtual vertex_list_graph_tag { };

    typedef CGAL::Triangulation_2<GT,TDS> Triangulation;

    typedef typename CGAL::Triangulation_2<GT,TDS>::Vertex_handle vertex_descriptor;
    typedef typename CGAL::Triangulation_2<GT,TDS>::Face_handle face_descriptor;
    typedef CGAL::detail::Edge<CGAL::Triangulation_2<GT,TDS>, typename CGAL::Triangulation_2<GT,TDS>::Edge>  edge_descriptor;
    typedef typename CGAL::Triangulation_2<GT,TDS>::All_edges_iterator  edge_iterator;


    typedef CGAL::detail::T2_halfedge_descriptor<Triangulation> halfedge_descriptor;

    typedef typename CGAL::Triangulation_2<GT,TDS>::All_halfedges_iterator  halfedge_iterator;

    typedef CGAL::detail::boost_all_vertices_iterator<Triangulation> vertex_iterator;
    typedef CGAL::detail::boost_all_faces_iterator<Triangulation> face_iterator;
    typedef CGAL::Counting_iterator<CGAL::detail::Out_edge_circulator<typename Triangulation::Edge_circulator, edge_descriptor>, edge_descriptor > out_edge_iterator;
    typedef CGAL::Counting_iterator<CGAL::detail::In_edge_circulator<typename Triangulation::Edge_circulator, edge_descriptor>, edge_descriptor > in_edge_iterator;
    typedef CGAL::Counting_iterator<typename Triangulation::Vertex_circulator> Incident_vertices_iterator;
    typedef Incident_vertices_iterator adjacency_iterator;

    typedef undirected_tag directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category; 
    typedef T2_graph_traversal_category traversal_category;
    typedef typename Triangulation::size_type size_type;
    typedef size_type vertices_size_type;
    typedef size_type edges_size_type;
    typedef size_type halfedges_size_type;
    typedef size_type faces_size_type;
    typedef size_type degree_size_type;

  // nulls
  static vertex_descriptor   null_vertex() { return vertex_descriptor(); }
  static face_descriptor     null_face()   { return face_descriptor(); }
  static halfedge_descriptor     null_halfedge()   { return halfedge_descriptor(); }
  };


} // namespace boost


namespace CGAL {

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor
  next(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor e,
       const Triangulation_2<Gt,Tds>& g)
  {
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor halfedge_descriptor;
    return halfedge_descriptor(e.first, g.ccw(e.second));
  }


  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor
  prev(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor e,
       const Triangulation_2<Gt,Tds>& g)
  {
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor halfedge_descriptor;
    return halfedge_descriptor(e.first, g.cw(e.second));
  }

  

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor
  opposite(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor e,
       const Triangulation_2<Gt,Tds>& g)
  {
    
    return g.mirror_edge(e);
  }
  
  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_descriptor
  source(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::edge_descriptor e,
         const Triangulation_2<Gt,Tds>& g)
  {
    return e.first->vertex(g.ccw(e.second));
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_descriptor
  target(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::edge_descriptor e,
         const Triangulation_2<Gt,Tds>& g)
  {
    return e.first->vertex(g.cw(e.second));
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_descriptor
  source(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor e,
         const Triangulation_2<Gt,Tds>& g)
  {
    return e.first->vertex(g.ccw(e.second));
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_descriptor
  target(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor e,
         const Triangulation_2<Gt,Tds>& g)
  {
    return e.first->vertex(g.cw(e.second));
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::face_descriptor
  face(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor e,
         const Triangulation_2<Gt,Tds>&)
  {
    return e.first;
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor
  halfedge(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::face_descriptor f,
           const Triangulation_2<Gt,Tds>&)
  {
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor halfedge_descriptor;
    return halfedge_descriptor(f,0);
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor
  halfedge(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_descriptor v,
           const Triangulation_2<Gt,Tds>& g)
  {
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >::face_descriptor face_descriptor;
    face_descriptor fd = v->face();
    int i = fd->index(v);
    return halfedge_descriptor(fd,g.ccw(i));
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor
  halfedge(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::edge_descriptor e,
           const Triangulation_2<Gt,Tds>&)
  {
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor halfedge_descriptor;
    return halfedge_descriptor(e.first,e.second);
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::edge_descriptor
  edge(typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_descriptor e,
           const Triangulation_2<Gt,Tds>&)
  {
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >::edge_descriptor edge_descriptor;
    return edge_descriptor(e.first,e.second);
  }


  template <class Gt, class Tds>
  inline Iterator_range<typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_iterator>  
  vertices(const Triangulation_2<Gt,Tds>& g)
  {
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_iterator
      Iter;
    return make_range( Iter(g.all_vertices_begin()), Iter(g.all_vertices_end()) );
  }


  template <class Gt, class Tds>
  inline Iterator_range<typename boost::graph_traits< Triangulation_2<Gt,Tds> >::edge_iterator>  
  edges(const Triangulation_2<Gt,Tds>& g)
  {    
    return make_range(g.all_edges_begin(), g.all_edges_end());
  }

  template <class Gt, class Tds>
  inline Iterator_range<typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedge_iterator >  
  halfedges(const Triangulation_2<Gt,Tds>& g)
  {    
    return make_range(g.all_halfedges_begin(), g.all_halfedges_end());
  }

  template <class Gt, class Tds>
  inline Iterator_range<typename boost::graph_traits< Triangulation_2<Gt,Tds> >::face_iterator >  
  faces(const Triangulation_2<Gt,Tds>& g)
  {
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >::face_iterator
      Iter;
    return make_range( Iter(g.all_faces_begin()), Iter(g.all_faces_end()) );
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::degree_size_type
  out_degree(
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const Triangulation_2<Gt,Tds>& g)
  {
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::degree_size_type deg = 0;
    typename Triangulation_2<Gt,Tds>::Edge_circulator c = g.incident_edges(u), done(c);
    if ( c != 0) {
        do {
            ++deg;
        } while (++c != done);
    }
    return deg;
  }

  template <class Gt, class Tds>
  inline Iterator_range<typename boost::graph_traits< Triangulation_2<Gt,Tds> >::out_edge_iterator >  
  out_edges(
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const Triangulation_2<Gt,Tds>& g)
  {
    typename Triangulation_2<Gt,Tds>::Edge_circulator ec(u,u->face());
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >
      ::out_edge_iterator Iter;
    
    return make_range( Iter(ec), Iter(ec,out_deg) );
  }

  template <class Gt, class Tds>
  inline Iterator_range<typename boost::graph_traits< Triangulation_2<Gt,Tds> >::in_edge_iterator >  
  in_edges(
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const Triangulation_2<Gt,Tds>& g)
  {
    typename Triangulation_2<Gt,Tds>::Edge_circulator ec(u,u->face());
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >
      ::in_edge_iterator Iter;
    return make_range( Iter(ec), Iter(ec,out_deg) );
  }

  template <class Gt, class Tds>
  inline Iterator_range<typename boost::graph_traits< Triangulation_2<Gt,Tds> >::adjacency_iterator>  
  adjacent_vertices(
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const Triangulation_2<Gt,Tds>& g)
  {
    typename Triangulation_2<Gt,Tds>::Vertex_circulator vc = out_edge_iterator(u,u.face());
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< Triangulation_2<Gt,Tds> >
      ::adjacency_iterator Iter;
    return make_range( Iter(vc), Iter(vc,out_deg) );
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertices_size_type
  num_vertices(const Triangulation_2<Gt,Tds>& g)
  {
    return g.tds().number_of_vertices();
  }  

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::edges_size_type
  num_edges(const Triangulation_2<Gt,Tds>& g)
  {

    return  g.tds().number_of_vertices() + g.tds().number_of_faces() - 2;
  }  

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::halfedges_size_type
  num_halfedges(const Triangulation_2<Gt,Tds>& g)
  {
    return  num_edges(g) * 2;
  }  

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::faces_size_type
  num_faces(const Triangulation_2<Gt,Tds>& g)
  {
    return  g.tds().number_of_faces();
  } 

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::degree_size_type
  in_degree(
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const Triangulation_2<Gt,Tds>& g)
  {
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::degree_size_type deg = 0;
    typename Triangulation_2<Gt,Tds>::Edge_circulator c = g.incident_edges(u), done(c);
    if ( c != 0) {
        do {
            ++deg;
        } while (++c != done);
    }
    return deg;
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< Triangulation_2<Gt,Tds> >::degree_size_type
  degree(
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const Triangulation_2<Gt,Tds>& g)
  {
    typename boost::graph_traits< Triangulation_2<Gt,Tds> >::degree_size_type deg = 0;
    typename Triangulation_2<Gt,Tds>::Edge_circulator c = g.incident_edges(u), done(c);
    if ( c != 0) {
        do {
            ++deg;
        } while (++c != done);
    }
    return deg;
  }


  // property maps
  template <class Gt, class Tds>
  class T2_vertex_id_map
    : public boost::put_get_helper<int, T2_vertex_id_map<Gt,Tds> >
  {
  public:
    typedef boost::readable_property_map_tag category;
    typedef int value_type;
    typedef int reference;
    typedef typename CGAL::Triangulation_2<Gt,Tds>::Vertex_handle key_type;
    
    T2_vertex_id_map()
    {}
    
    long operator[](key_type vh) const {
      return vh->id(); 
    }
  };

  template <class Gt, class Tds>
  class T2_vertex_point_map
  {
  public:
    typedef boost::lvalue_property_map_tag category;
    typedef typename Tds::Vertex::Point value_type;
    typedef value_type& reference;
    typedef typename CGAL::Triangulation_2<Gt,Tds>::Vertex_handle key_type;

    friend reference get(T2_vertex_point_map<Gt,Tds>, key_type vh)
    { 
      return vh->point(); 
    }
    friend void put(T2_vertex_point_map<Gt,Tds>, key_type vh, reference v)
    {
      vh->point()=v; 
    }
    reference operator[](key_type vh) const {
      return vh->point();
    }
  };


  template <class Gt, class Tds>
  class T2_edge_id_map
    : public boost::put_get_helper<int, T2_edge_id_map<Gt,Tds> >
  {
  public:
    typedef boost::readable_property_map_tag category;
    typedef int value_type;
    typedef int reference;
    typedef typename CGAL::Triangulation_2<Gt,Tds>::Edge key_type;
    
    T2_edge_id_map()
    {}
    
    long operator[](key_type e) const {
      return (3 * e.first.id()) + e.second; 
    }
  };

  template <class Gt, class Tds>
  class T2_edge_weight_map
    : public boost::put_get_helper<typename Gt::FT, T2_edge_weight_map<Gt, Tds> >
  {
  private:
    const CGAL::Triangulation_2<Gt,Tds>& tr;
  public:
    typedef boost::readable_property_map_tag category;
    typedef typename Gt::FT value_type;
    typedef value_type reference;
    typedef typename CGAL::Triangulation_2<Gt,Tds>::Edge key_type;

    T2_edge_weight_map(const CGAL::Triangulation_2<Gt,Tds>& tr_) 
      : tr(tr_) 
    { }

    value_type operator[](key_type e) const {
      return tr.segment(e).squared_length();
    }
  };


  template <class Gt, class Tds>
  inline T2_vertex_id_map<Gt,Tds>
  get(boost::vertex_index_t, const Triangulation_2<Gt,Tds>&) {
    T2_vertex_id_map<Gt,Tds> m;
    return m;
  }

  template <class Gt, class Tds>
  inline T2_vertex_point_map<Gt,Tds>
  get(boost::vertex_point_t, const Triangulation_2<Gt,Tds>&) {
    T2_vertex_point_map<Gt,Tds> m;
    return m;
  }

  template <class Gt, class Tds>
  inline T2_edge_id_map<Gt,Tds>
  get(boost::edge_index_t, const Triangulation_2<Gt,Tds>&) {
    T2_edge_id_map<Gt,Tds> m;
    return m;
  }

  template <class Gt, class Tds>
  inline T2_edge_weight_map<Gt,Tds>
  get(boost::edge_weight_t, const Triangulation_2<Gt,Tds>& g) {
    T2_edge_weight_map<Gt,Tds> m(g);
    return m;
  }

  template <class Tag>
  struct T2_property_map { };

  template <>
  struct T2_property_map<boost::vertex_index_t> {
    template <class Gt, class Tds>
    struct bind_ {
      typedef T2_vertex_id_map<Gt,Tds> type;
      typedef T2_vertex_id_map<Gt,Tds> const_type;
    };
  };



  template <>
  struct T2_property_map<boost::vertex_point_t> {
    template <class Gt, class Tds>
    struct bind_ {
      typedef T2_vertex_point_map<Gt,Tds> type;
      typedef T2_vertex_point_map<Gt,Tds> const_type;
    };
  };


  template <>
  struct T2_property_map<boost::edge_index_t> {
    template <class Gt, class Tds>
    struct bind_ {
      typedef T2_edge_id_map<Gt,Tds> type;
      typedef T2_edge_id_map<Gt,Tds> const_type;
    };
  };


  template <>
  struct T2_property_map<boost::edge_weight_t> {
    template <class Gt, class Tds>
    struct bind_ {
      typedef T2_edge_weight_map<Gt,Tds> type;
      typedef T2_edge_weight_map<Gt,Tds> const_type;
    };
  };

} // namespace CGAL

namespace boost {
  // g++ 'enumeral_type' in template unification not implemented workaround
  template <class Gt, class Tds, class Tag>
  struct property_map<CGAL::Triangulation_2<Gt,Tds>, Tag> {
    typedef typename 
    CGAL::T2_property_map<Tag>::template bind_<Gt,Tds> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };

  // see struct property_map in Polyehdron for an explanation
  template <class Gt, class Tds, class Tag>
  struct property_map<const CGAL::Triangulation_2<Gt,Tds>, Tag> {
    typedef typename 
    CGAL::T2_property_map<Tag>::template bind_<Gt,Tds> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };

} // namespace boost


namespace CGAL {

  template <class Gt, class Tds, class PropertyTag, class Key>
  inline
  typename boost::property_traits<
    typename boost::property_map<Triangulation_2<Gt,Tds>,PropertyTag>::const_type>::value_type
  get(PropertyTag p, const Triangulation_2<Gt,Tds>& g, const Key& key) {
    return get(get(p, g), key);
  }
  
  template <class Gt, class Tds, class PropertyTag, class Key,class Value>
  inline void
  put(PropertyTag p, Triangulation_2<Gt,Tds>& g, 
      const Key& key, const Value& value)
  {
    typedef typename boost::property_map<Triangulation_2<Gt,Tds>, PropertyTag>::type Map;
    Map pmap = get(p, g);
    put(pmap, key, value);
  }

} // namespace CGAL 

namespace boost {

  // What are those needed for ???
  template <typename Gt, typename Tds>
  struct edge_property_type<CGAL::Triangulation_2<Gt,Tds> > {
    typedef void type;
  };  

  template <typename Gt, typename Tds>
  struct vertex_property_type<CGAL::Triangulation_2<Gt,Tds> > {
    typedef void type;
  };
} // namespace boost


namespace std {


#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4099) // For VC10 it is class hash 
#endif

#ifndef CGAL_CFG_NO_STD_HASH

  template < class T, class EdgeBase>
  struct hash<CGAL::detail::Edge<T,EdgeBase> > {
    std::size_t operator()(const CGAL::detail::Edge<T,EdgeBase>& e) const
    {
      std::cerr << "Triangulation_2::Edge HashFct" << std::endl;
      return hash_value(e);
    }
  }; 

  template < class Tr>
  struct hash<CGAL::detail::T2_halfedge_descriptor<Tr> > {
    std::size_t operator()(const CGAL::detail::T2_halfedge_descriptor<Tr>& e) const
    {
      std::cerr << "Triangulation_2::halfedge_descriptor HashFct" << std::endl;
      return hash_value(e);
    }
  };

#endif // CGAL_CFG_NO_STD_HASH

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // namespace std

//#include <CGAL/graph_traits_Delaunay_triangulation_2.h>

#endif // CGAL_GRAPH_TRAITS_TRIANGULATION_2_H
