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

#ifndef CGAL_POLYHEDRON_EXTENDED_BGL_H
#define CGAL_POLYHEDRON_EXTENDED_BGL_H

#include <CGAL/HalfedgeDS_items_decorator.h>
#include <CGAL/Polyhedron_BGL.h>

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
#  define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS
#else
#  define CGAL_HDS_PARAM_ class HDS
#endif

CGAL_BEGIN_NAMESPACE

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edges_size_type 
num_undirected_edges(const Polyhedron_3<Gt,I,HDS,A>& p)
{
  return p.size_of_halfedges() / 2 ;
}

//
// Const versions
// 

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline std::pair<typename boost::undirected_graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_iterator
                ,typename boost::undirected_graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_iterator 
                >  
undirected_edges( Polyhedron_3<Gt,I,HDS,A> const& p )
{
  typedef typename boost::undirected_graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_iterator Iter;
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

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor
out_edge( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor u
        , Polyhedron_3<Gt,I,HDS,A> const& 
        )
{
  HalfedgeDS_items_decorator< Polyhedron_3<Gt,I,HDS,A> > D ;
  return D.get_vertex_edge(u)->opposite();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor
in_edge( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor u
       , Polyhedron_3<Gt,I,HDS,A> const& 
       )
{
  HalfedgeDS_items_decorator< Polyhedron_3<Gt,I,HDS,A> > D ;
  return D.get_vertex_edge(u);
}

//
// Non-Const versions
// 

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline std::pair<typename boost::undirected_graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_iterator
                ,typename boost::undirected_graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_iterator 
                >  
undirected_edges( Polyhedron_3<Gt,I,HDS,A>& p )
{
  typedef typename boost::undirected_graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_iterator Iter;
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



template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor 
out_edge( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor u
        , Polyhedron_3<Gt,I,HDS,A>& 
        )
{
  HalfedgeDS_items_decorator< Polyhedron_3<Gt,I,HDS,A> > D ;
  return D.get_vertex_edge(u)->opposite();
}


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor 
in_edge( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor u
       , Polyhedron_3<Gt,I,HDS,A>& 
       )
{
  HalfedgeDS_items_decorator< Polyhedron_3<Gt,I,HDS,A> > D ;
  return D.get_vertex_edge(u);
}

CGAL_END_NAMESPACE

#undef CGAL_HDS_

#endif // CGAL_POLYHEDRON_EXTENDED_BGL_H
