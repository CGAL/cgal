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

#ifndef CGAL_BOOST_GRAPH_GRAPH_TRAITS_HALFEDGEDS_H
#define CGAL_BOOST_GRAPH_GRAPH_TRAITS_HALFEDGEDS_H

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/basic.h>
#include <CGAL/Counting_iterator.h>

namespace CGAL {

template <class Circ, class E>
class HDS_in_halfedge_circulator : public Circ
{
private:
  mutable E e;

public:

  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  HDS_in_halfedge_circulator()
    : Circ()
  {}

  HDS_in_halfedge_circulator(Circ c)
    : Circ(c)
  {}

  const E& operator*() const
  {
    e = *this;
    return e;
  }
};

template <class Circ, class E>
class HDS_out_halfedge_circulator : public Circ
{
private:
  mutable E e;

public:

  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  HDS_out_halfedge_circulator()
    : Circ()
  {}

  HDS_out_halfedge_circulator(Circ c)
    : Circ(c)
  {}

  const E& operator*() const
  {
    e = *this;
    e = e->opposite();
    return e;
  }
};
  

//  The vertex iterator of the bgl must evaluate to a vertex handle, not to a vertex
template < class HDS, class Vertex_iterator, class Vertex_handle>
class HDS_all_vertices_iterator_base {
protected:
    Vertex_iterator nt;
public:
  typedef Vertex_iterator  Iterator;
  typedef HDS_all_vertices_iterator_base<HDS,Vertex_iterator,Vertex_handle> Self;

  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;
  typedef typename std::iterator_traits<Iterator>::difference_type   difference_type;
  typedef Vertex_handle                                              value_type;
  typedef value_type                                                 reference;
  typedef value_type                                                 pointer;

protected:

  HDS_all_vertices_iterator_base() {}
  HDS_all_vertices_iterator_base( Iterator j) : nt(j) {}

public:

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
};

template < class HDS >
class HDS_all_vertices_const_iterator 
  : public HDS_all_vertices_iterator_base<HDS,typename HDS::Vertex_const_iterator,typename HDS::Vertex_const_handle>
{
  typedef HDS_all_vertices_iterator_base<HDS,typename HDS::Vertex_const_iterator,typename HDS::Vertex_const_handle> Base ;
  
public:

  typedef typename HDS::Vertex_const_iterator Iterator;

  HDS_all_vertices_const_iterator() {}
  HDS_all_vertices_const_iterator( Iterator j) : Base(j) {}
};


template < class HDS >
class HDS_all_vertices_iterator 
  : public HDS_all_vertices_iterator_base<HDS,typename HDS::Vertex_iterator,typename HDS::Vertex_handle>
{
  typedef HDS_all_vertices_iterator_base<HDS,typename HDS::Vertex_iterator,typename HDS::Vertex_handle> Base ;
  
public:

  typedef typename HDS::Vertex_iterator Iterator;

  HDS_all_vertices_iterator() {}
  HDS_all_vertices_iterator( Iterator j) : Base(j) {}
};

template < class HDS, class Iterator_, class Value_type>
class HDS_all_edges_iterator_base {
protected:
  Iterator_ nt;
public:
  typedef Iterator_  Iterator;
  typedef HDS_all_edges_iterator_base<HDS,Iterator_,Value_type> Self;

  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;
  typedef typename std::iterator_traits<Iterator>::difference_type   difference_type;
  typedef Value_type                                                 value_type;
  typedef value_type                                                 reference;
  typedef value_type                                                 pointer;
  
protected:

  HDS_all_edges_iterator_base() {}
  HDS_all_edges_iterator_base( Iterator j) : nt(j) {}

public:

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
};

template < class HDS >
class HDS_all_halfedges_const_iterator 
  : public HDS_all_edges_iterator_base<HDS,typename HDS::Halfedge_const_iterator,typename HDS::Halfedge_const_handle>
{
  typedef HDS_all_edges_iterator_base<HDS,typename HDS::Halfedge_const_iterator,typename HDS::Halfedge_const_handle> Base ;
  
public:

  typedef typename HDS::Halfedge_const_iterator Iterator;

  HDS_all_halfedges_const_iterator() {}
  HDS_all_halfedges_const_iterator( Iterator j) : Base(j) {}
};

template < class HDS >
class HDS_all_halfedges_iterator 
  : public HDS_all_edges_iterator_base<HDS,typename HDS::Halfedge_iterator,typename HDS::Halfedge_handle>
{
  typedef HDS_all_edges_iterator_base<HDS,typename HDS::Halfedge_iterator,typename HDS::Halfedge_handle> Base ;
  
public:

  typedef typename HDS::Halfedge_iterator Iterator;

  HDS_all_halfedges_iterator() {}
  HDS_all_halfedges_iterator( Iterator j) : Base(j) {}
};

template <class HDS_>
struct HDS_graph_traits
{
public :
  
  struct HDS_graph_traversal_category : public virtual boost::bidirectional_graph_tag,
                                        public virtual boost::vertex_list_graph_tag,
                                        public virtual boost::edge_list_graph_tag
  {};

  typedef HDS_ HDS;

  typedef typename HDS::Vertex_handle   vertex_descriptor;
  typedef typename HDS::Halfedge_handle edge_descriptor;
  
  typedef HDS_all_vertices_iterator<HDS>  vertex_iterator;
  typedef HDS_all_halfedges_iterator<HDS>  edge_iterator;
  
private:

  typedef typename HDS::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator ;
  
  typedef HDS_out_halfedge_circulator<Halfedge_around_vertex_circulator,edge_descriptor> out_edge_circulator ;
  typedef HDS_in_halfedge_circulator <Halfedge_around_vertex_circulator,edge_descriptor> in_edge_circulator ;
  
public :

  typedef Counting_iterator<out_edge_circulator, edge_descriptor> out_edge_iterator;
  typedef Counting_iterator<in_edge_circulator , edge_descriptor> in_edge_iterator;
                                 
  typedef boost::directed_tag               directed_category;
  typedef boost::disallow_parallel_edge_tag edge_parallel_category; 
  typedef HDS_graph_traversal_category      traversal_category;
  
  typedef typename HDS::size_type vertices_size_type;
  typedef vertices_size_type      edges_size_type;
  typedef vertices_size_type      degree_size_type;
};


} //namespace CGAL

#endif // CGAL_BOOST_GRAPH_GRAPH_TRAITS_HALFEDGEDS_H
