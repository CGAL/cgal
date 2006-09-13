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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Polyhedron_extended_BGL.h $
// $Id: Polyhedron_extended_BGL.h 32155 2006-06-29 17:08:41Z fcacciola $
// 
//
// Author(s) : Fernando Caccciola <fernando.cacciola@gmail.com>

#ifndef CGAL_BOOST_GRAPH_POLYHEDRON_HALFEDGE_GRAPH_TRAITS_H
#define CGAL_BOOST_GRAPH_POLYHEDRON_HALFEDGE_GRAPH_TRAITS_H

#include <CGAL/HalfedgeDS_items_decorator.h>
#include <CGAL/boost/graph/Polyhedron_BGL.h>
#include <CGAL/boost/graph/Extended_BGL.h>

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
#  define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS
#else
#  define CGAL_HDS_PARAM_ class HDS
#endif

CGAL_BEGIN_NAMESPACE


//
// Const versions
// 
template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct halfedge_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const > 
  : CGAL::HDS_halfedge_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>
{};


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline std::pair<typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::undirected_edge_iterator
                ,typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::undirected_edge_iterator 
                >  
undirected_edges( Polyhedron_3<Gt,I,HDS,A> const& p )
{
  typedef typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::undirected_edge_iterator Iter;
  return std::make_pair( Iter(p.edges_begin()), Iter(p.edges_end()) );
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor
next_edge( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor outedge
         , Polyhedron_3<Gt,I,HDS,A> const& 
         )
{
  return outedge->next();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor 
prev_edge( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor outedge
         , Polyhedron_3<Gt,I,HDS,A> const& 
         )
{
  HalfedgeDS_items_decorator< Polyhedron_3<Gt,I,HDS,A> > D ;
  return D.get_prev(outedge);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor
opposite_edge( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor e
             , Polyhedron_3<Gt,I,HDS,A> const& 
             )
{
  return e->opposite();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor
next_edge_ccw( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor outedge
             , Polyhedron_3<Gt,I,HDS,A> const& 
             )
{
  HalfedgeDS_items_decorator< Polyhedron_3<Gt,I,HDS,A> > D ;
  return D.get_prev(outedge)->opposite();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor
next_edge_cw( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor outedge
            , Polyhedron_3<Gt,I,HDS,A> const& 
            )
{
  return outedge->opposite()->next();
}

//
// Non-Const versions
// 
template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct geometric_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> > 
{
  typedef CGAL::Polyhedron_3<Gt,I,HDS,A>  Polyhedron ;
  
  typedef typename Polyhedron::Point_3 Point ;
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct halfedge_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> > 
  : CGAL::HDS_halfedge_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >
{};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline std::pair<typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> >::undirected_edge_iterator
                ,typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> >::undirected_edge_iterator 
                >  
undirected_edges( Polyhedron_3<Gt,I,HDS,A>& p )
{
  typedef typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> >::undirected_edge_iterator Iter;
  return std::make_pair( Iter(p.edges_begin()), Iter(p.edges_end()) );
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor 
next_edge( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor outedge
         , Polyhedron_3<Gt,I,HDS,A>& 
         )
{
  return outedge->next();
}


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor 
prev_edge( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor outedge
         , Polyhedron_3<Gt,I,HDS,A>& p
         )
{
  CGAL::HalfedgeDS_items_decorator< Polyhedron_3<Gt,I,HDS,A> > D ;
  return D.get_prev(outedge);
}


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor
opposite_edge( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor e
             , Polyhedron_3<Gt,I,HDS,A>& 
             )
{
  return e->opposite();
}


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor 
next_edge_ccw( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor outedge
             , Polyhedron_3<Gt,I,HDS,A>& 
             )
{
  HalfedgeDS_items_decorator< Polyhedron_3<Gt,I,HDS,A> > D ;
  return D.get_prev(outedge)->opposite();
}


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor
next_edge_cw( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor outedge
            , Polyhedron_3<Gt,I,HDS,A>& 
            )
{
  return outedge->opposite()->next();
}


CGAL_END_NAMESPACE

#undef CGAL_HDS_

#endif // CGAL_BOOST_GRAPH_POLYHEDRON_HALFEDGE_GRAPH_TRAITS_H
