// Copyright (c) 2006 Geometry Factory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Andreas Fabri  <andreas.fabri@geometryfactory.com>)

#ifndef CGAL_POLYHEDRON_GRAPH_TRAITS_3_H
#define CGAL_POLYHEDRON_GRAPH_TRAITS_3_H

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/HalfedgeDS_items_decorator.h>


namespace boost { 

  namespace detail {



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
  e = *this;
    return e;
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
    e = *this;
    e = e->opposite();
    return e;
  }
};
  

//  The vertex iterator of the bgl must evaluate to a vertex handle, not to a vertex
template < class P, class Vertex_iterator, class Vertex_handle>
class boost_P_all_vertices_iterator_base {
protected:
    Vertex_iterator nt;
public:
  typedef Vertex_iterator  Iterator;
  typedef boost_P_all_vertices_iterator_base<P,Vertex_iterator,Vertex_handle> Self;

  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;
  typedef Vertex_handle  value_type;
  typedef typename std::iterator_traits<Iterator>::difference_type           difference_type;
  typedef value_type      reference;
  typedef value_type      pointer;

protected:

  // CREATION
  // --------

  boost_P_all_vertices_iterator_base()
  {}

  boost_P_all_vertices_iterator_base( Iterator j) : nt(j) {}

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

template < class P >
class boost_P_all_vertices_const_iterator 
  : public boost_P_all_vertices_iterator_base<P,typename P::Vertex_const_iterator,typename P::Vertex_const_handle>
{
  typedef boost_P_all_vertices_iterator_base<P,typename P::Vertex_const_iterator,typename P::Vertex_const_handle> Base ;
  
public:

  typedef typename P::Vertex_const_iterator Iterator;

  // CREATION
  // --------

  boost_P_all_vertices_const_iterator() {}

  boost_P_all_vertices_const_iterator( Iterator j) : Base(j) {}
};

template < class P >
class boost_P_all_vertices_iterator 
  : public boost_P_all_vertices_iterator_base<P,typename P::Vertex_iterator,typename P::Vertex_handle>
{
  typedef boost_P_all_vertices_iterator_base<P,typename P::Vertex_iterator,typename P::Vertex_handle> Base ;
  
public:

  typedef typename P::Vertex_iterator Iterator;

  // CREATION
  // --------

  boost_P_all_vertices_iterator() {}

  boost_P_all_vertices_iterator( Iterator j) : Base(j) {}
};

template < class P, class Halfedge_iterator, class Halfedge_handle >
class boost_P_all_edges_iterator_base {
protected:
  Halfedge_iterator nt;
public:
  typedef Halfedge_iterator  Iterator;
  typedef boost_P_all_edges_iterator_base<P,Halfedge_iterator,Halfedge_handle> Self;

  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;
  typedef Halfedge_handle  value_type;
  typedef typename std::iterator_traits<Iterator>::difference_type           difference_type;
  typedef value_type      reference;
  typedef value_type      pointer;
  
protected:

  // CREATION
  // --------

  boost_P_all_edges_iterator_base()
  {}

  boost_P_all_edges_iterator_base( Iterator j) : nt(j) {}

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

template < class P >
class boost_P_all_edges_const_iterator 
  : public boost_P_all_edges_iterator_base<P,typename P::Edge_const_iterator,typename P::Halfedge_const_handle>
{
  typedef boost_P_all_edges_iterator_base<P,typename P::Edge_const_iterator,typename P::Halfedge_const_handle> Base ;
  
public:

  typedef typename P::Edge_const_iterator Iterator;

  // CREATION
  // --------

  boost_P_all_edges_const_iterator() {}

  boost_P_all_edges_const_iterator( Iterator j) : Base(j) {}
};

template < class P >
class boost_P_all_edges_iterator 
  : public boost_P_all_edges_iterator_base<P,typename P::Edge_iterator,typename P::Halfedge_handle>
{
  typedef boost_P_all_edges_iterator_base<P,typename P::Edge_iterator,typename P::Halfedge_handle> Base ;
  
public:

  typedef typename P::Edge_iterator Iterator;

  // CREATION
  // --------

  boost_P_all_edges_iterator() {}

  boost_P_all_edges_iterator( Iterator j) : Base(j) {}
};


  } // namespace detail 

template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  struct graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> > {

    struct Polyhedron_graph_traversal_category : 
      public virtual bidirectional_graph_tag,
      public virtual adjacency_graph_tag,
      public virtual vertex_list_graph_tag,
      public virtual edge_list_graph_tag
       { };

    typedef CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> Polyhedron;

    typedef typename Polyhedron::Vertex_handle   vertex_descriptor;
    typedef typename Polyhedron::Halfedge_handle edge_descriptor;

    // NOTE: These are NOT defined as part of the BGL graph_traits class
    typedef typename Polyhedron::Vertex_const_handle   vertex_const_descriptor;
    typedef typename Polyhedron::Halfedge_const_handle edge_const_descriptor;
    
    typedef detail::boost_P_all_vertices_iterator<Polyhedron>  vertex_iterator;
    typedef detail::boost_P_all_edges_iterator   <Polyhedron>  edge_iterator;
    
    // NOTE: These are NOT defined as part of the BGL graph_traits class
    typedef detail::boost_P_all_vertices_const_iterator<Polyhedron>  vertex_const_iterator;
    typedef detail::boost_P_all_edges_const_iterator   <Polyhedron>  edge_const_iterator;
    
    typedef CGAL::Counting_iterator<detail::Out_edge_circulator<typename Polyhedron::Halfedge_around_vertex_circulator
                                                               ,edge_descriptor
                                                               >
                                   , edge_descriptor
                                   > out_edge_iterator;
                                   
    typedef CGAL::Counting_iterator<detail::In_edge_circulator<typename Polyhedron::Halfedge_around_vertex_circulator
                                                              ,edge_descriptor
                                                              >
                                   , edge_descriptor
                                   > in_edge_iterator;
                                   
    // NOTE: This is NOT defined as part of the BGL graph_traits class
    typedef CGAL::Counting_iterator<detail::Out_edge_circulator<typename Polyhedron::Halfedge_around_vertex_const_circulator
                                                               ,edge_const_descriptor
                                                               >
                                   , edge_const_descriptor
                                   > out_edge_const_iterator;
                                   
    // NOTE: This is NOT defined as part of the BGL graph_traits class
    typedef CGAL::Counting_iterator<detail::In_edge_circulator<typename Polyhedron::Halfedge_around_vertex_const_circulator
                                                              ,edge_const_descriptor
                                                              >
                                   , edge_const_descriptor
                                   > in_edge_const_iterator;
                                   
    //typedef CGAL::Counting_iterator<typename Triangulation::Vertex_circulator> Incident_vertices_iterator;
    //typedef Incident_vertices_iterator adjacency_iterator;

    typedef directed_tag directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category; 
    typedef Polyhedron_graph_traversal_category traversal_category;
    
    typedef typename Polyhedron::size_type vertices_size_type;
    typedef typename Polyhedron::size_type edges_size_type;
    typedef typename Polyhedron::size_type degree_size_type;
  };


  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_descriptor
  source(typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_descriptor e,
         CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> & p)
  {
    return e->opposite()->vertex();
  }
  
  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_const_descriptor
  source(typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_descriptor e,
         CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const& p)
  {
    return e->opposite()->vertex();
  }


  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_descriptor
  target(typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_descriptor e,
         CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> & p)
  {
    return e->vertex();
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_const_descriptor
  target(typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_descriptor e,
         CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const& p)
  {
    return e->vertex();
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  inline std::pair<
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_iterator,
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_iterator >  
  vertices( CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& p)
  {
    typedef typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_iterator
      Iter;
    return std::make_pair( Iter(p.vertices_begin()), Iter(p.vertices_end()) );
  }
  
  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  inline std::pair<
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_const_iterator,
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_const_iterator >  
  vertices( CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const& p)
  {
    typedef typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_const_iterator
      Iter;
    return std::make_pair( Iter(p.vertices_begin()), Iter(p.vertices_end()) );
  }


  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  inline std::pair<
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_iterator,
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_iterator >  
  edges( CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& p)
  {
    typedef typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_iterator
      Iter;
    return std::make_pair( Iter(p.edges_begin()), Iter(p.edges_end()) );
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  inline std::pair<
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_iterator,
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_iterator >  
  edges( CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const& p)
  {
    typedef typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_iterator
      Iter;
    return std::make_pair( Iter(p.edges_begin()), Iter(p.edges_end()) );
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  inline std::pair<
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::in_edge_iterator,
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::in_edge_iterator >  
  in_edges(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_descriptor u, 
    CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& g)
  {
    typename CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>::Halfedge_around_vertex_circulator 
      ec = u->vertex_begin();
    int in_deg = in_degree(u,g);
    typedef typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >
      ::in_edge_iterator Iter;
    return std::make_pair( Iter(ec), Iter(ec,in_deg) );
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  inline std::pair<
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::in_edge_const_iterator,
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::in_edge_const_iterator >  
  in_edges(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_const_descriptor u, 
    CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const& g)
  {
    typename CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>::Halfedge_around_vertex_const_circulator 
      ec = u->vertex_begin();
    int in_deg = in_degree(u,g);
    typedef typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >
      ::in_edge_const_iterator Iter;
    return std::make_pair( Iter(ec), Iter(ec,in_deg) );
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  inline std::pair<
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::out_edge_iterator,
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::out_edge_iterator >  
  out_edges(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_descriptor u, 
    const CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& g)
  {
    typename CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>::Halfedge_around_vertex_circulator
       ec = u->vertex_begin();
    int out_deg = out_degree(u,g);
    typedef typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >
      ::out_edge_iterator Iter;
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  inline std::pair<
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::out_edge_const_iterator,
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::out_edge_const_iterator >  
  out_edges(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_const_descriptor u, 
    CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const& g)
  {
    typename CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>::Halfedge_around_vertex_const_circulator
       ec = u->vertex_begin();
    int out_deg = out_degree(u,g);
    typedef typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >
      ::out_edge_const_iterator Iter;
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertices_size_type
  num_vertices(const CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& p)
  {
    return p.size_of_vertices();
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edges_size_type
  num_edges(const CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& p)
  {
    return p.size_of_halfedges() / 2 ;
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::degree_size_type
  degree(typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_const_descriptor v,
         const CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>&)
  {
    return v->vertex_degree();
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::degree_size_type
  out_degree(typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_const_descriptor v,
         const CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>&)
  {
    return v->vertex_degree();
  }

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::degree_size_type
  in_degree(typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_const_descriptor v,
         const CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>&)
  {
    return v->vertex_degree();
  }


 template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  class Polyhedron_edge_weight_map
    : public put_get_helper<double, Polyhedron_edge_weight_map<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >
  {
  private:
  
    typedef CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> Polyhedron ;
    
    Polyhedron const* p;
    
  public:
  
    typedef readable_property_map_tag category;
    typedef double value_type;
    typedef double reference;
    typedef typename graph_traits<Polyhedron>::edge_const_descriptor  key_type;

    Polyhedron_edge_weight_map( Polyhedron const& p_) 
      : p( boost::addressof(p_) ) 
    { }

    reference operator[](key_type const& e) const {
      return CGAL::squared_distance(e->vertex()->point(), e->opposite()->vertex()->point());
    }
  };


  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  inline Polyhedron_edge_weight_map<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>
  get(edge_weight_t, const CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& p) {
    Polyhedron_edge_weight_map<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> m(p);
    return m;
  }

 template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  class Polyhedron_vertex_point_map
    : public put_get_helper< typename PolyhedronTraits_3::Point_3
                          , Polyhedron_vertex_point_map<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> 
                          >
  {
  private:
  
    typedef CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> Polyhedron ;
  
    Polyhedron* p;
    
  public:
  
    typedef typename PolyhedronTraits_3::Point_3 Point_3 ;
    
    typedef lvalue_property_map_tag category;
    
    typedef Point_3   value_type;
    typedef Point_3 & reference;
    
    typedef typename graph_traits<Polyhedron>::vertex_descriptor key_type;

    Polyhedron_vertex_point_map( Polyhedron& p_) 
      : p( boost::addressof(p_) ) 
    { }

    reference operator[](key_type const& e) const {
      return e->point();
    }
  };
  
  struct vertex_point_t {} ;
  
  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  inline Polyhedron_vertex_point_map<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>
  get(vertex_point_t, CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& p) {
    Polyhedron_vertex_point_map<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> m(p);
    return m;
  }
  
  template <class Tag>
  struct Polyhedron_property_map { };

  
  template <>
  struct Polyhedron_property_map<edge_weight_t> {
    template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
    struct bind_ {
      typedef Polyhedron_edge_weight_map<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> type;
      typedef Polyhedron_edge_weight_map<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const_type;
    };
  };

  
  template <>
  struct Polyhedron_property_map<vertex_point_t> {
    template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
    struct bind_ {
      typedef Polyhedron_vertex_point_map<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> type;
      typedef Polyhedron_vertex_point_map<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const_type;
    };
  };

  // g++ 'enumeral_type' in template unification not implemented workaround
  template <class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc, class Tag>
  struct property_map<CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>, Tag> {
    typedef typename 
      Polyhedron_property_map<Tag>::template bind_<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };

  template <class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc, class PropertyTag, class Key>
  inline
  typename boost::property_traits<
    typename boost::property_map<CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>,PropertyTag>::type>::value_type
  get(PropertyTag p, CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& g, const Key& key) {
    return get(get(p, g), key);
  }
  
  template <class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc, class PropertyTag, class Key>
  inline
  typename boost::property_traits<
    typename boost::property_map<CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>,PropertyTag>::const_type>::value_type
  get(PropertyTag p, const CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& g, const Key& key) {
    return get(get(p, g), key);
  }
  
  template <class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc, class PropertyTag, class Key,class Value>
  inline void
  put(PropertyTag p, CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& g, 
      const Key& key, const Value& value)
  {
    typedef typename property_map<CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>, PropertyTag>::type Map;
    Map pmap = get(p, g);
    put(pmap, key, value);
  }


  // What are those needed for ???
  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  struct edge_property_type<CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> > {
    typedef edge_weight_t type;
  };  

  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  struct vertex_property_type<CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> > {
    typedef vertex_point_t type;
  };


  //****************************************************************************************************************
  //
  //                                           Extended-BGL interface
  //  
  //
  // Some queries typical of CGAL data-structures can be handled by the BGL interface but with an unneccesary
  // high complexity. For example, accessing the halfedge incident upon a HDS vertex can be done via out_edges(), 
  // but out_edges() has complexity O(number_of_adjacent_vertices).
  // For that reason, the following "Extended-BGL" interface is defined.
  // Algorithms MUST indicate whether they need a DS to support this Extended-BGL interface.
  //
  // first_in_edge(v,dag)
  // first_out_edge(v,dag)
  // next_in_edge_ccw(e,dag)
  // next_in_edge_cw(e,dag)
  // next_out_edge_ccw(e,dag)
  // next_out_edge_cw(e,dag)
  // opposite_edge(e,dag) 
  //
  
  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_descriptor
  out_edge(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_descriptor u, 
    CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& g)
  {
    CGAL::HalfedgeDS_items_decorator< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> > D ;
    return D.get_vertex_halfedge(u)->opposite();
  }
  
  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_descriptor
  out_edge(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::vertex_const_descriptor u, 
    CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const& g)
  {
    CGAL::HalfedgeDS_items_decorator< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> > D ;
    return D.get_vertex_halfedge(u)->opposite();
  }


    template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_descriptor
  next_edge_ccw(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_descriptor outedge, 
    CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& g)
  {
    CGAL::HalfedgeDS_items_decorator< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> > D ;
    return D.get_prev(outedge)->opposite();
  }
  
  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_descriptor
  next_edge_ccw(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_descriptor outedge, 
    CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const& g)
  {
    CGAL::HalfedgeDS_items_decorator< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> > D ;
    return D.get_prev(outedge)->opposite();
  }
  
    template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_descriptor
  next_edge_cw(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_descriptor outedge, 
    CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& g)
  {
    return outedge->opposite()->next();
  }
  
  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_descriptor
  next_edge_cw(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_descriptor outedge, 
    CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const& g)
  {
    return outedge->opposite()->next();
  }

      template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_descriptor
  opposite_edge(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_descriptor e, 
    CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& g)
  {
    return e->opposite();
  }
  
  template < class PolyhedronTraits_3, class PolyhedronItems_3, 
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
class T_HDS, class Alloc>
  typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_descriptor
  opposite_edge(
    typename graph_traits< CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> >::edge_const_descriptor e, 
    CGAL::Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> const& g)
  {
    return e->opposite();
  }
      
} // namespace boost

#endif // CGAL_POLYHEDRON_GRAPH_TRAITS_3_H
