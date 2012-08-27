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


#ifndef CGAL_BOOST_GRAPH_HALFEDGE_GRAPH_TRAITS_POLYHEDRON_3_H
#define CGAL_BOOST_GRAPH_HALFEDGE_GRAPH_TRAITS_POLYHEDRON_3_H

#include <CGAL/HalfedgeDS_items_decorator.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits.h>
#include <CGAL/boost/graph/halfedge_graph_traits_HalfedgeDS.h>

#define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS

//
// NOTE: The BGL algorithms are NOT const-correct: i.e., they take a "G const&" 
// but instantiate "graph_traits<G>" instead of "graph_traits<G const>"
// This is known Boost bug which will eventually be fixed, but in the meantime we need
// to coerce both const and non-const specializations.
// That is, HDS_graph_traits<G const> is really the same as HDS_graph_traits<G>
// so graph_traits<G> is also the same as graph_traits<G const>.
// Therefore, while, for instance, "graph_traits<G const>::vertex_descriptor"
// is conceptually const-correct, it actually corresponds to the non-const handle,
// hence the const_cast<> used below in the functions implementation.
//

namespace CGAL {

//
// Const versions
// 
template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct halfedge_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const > 
  : CGAL::HDS_halfedge_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> > // See NOTE above!
{
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline std::pair<typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::undirected_edge_iterator
                ,typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::undirected_edge_iterator 
                >  
undirected_edges( Polyhedron_3<Gt,I,HDS,A> const& p )
{
  typedef typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::undirected_edge_iterator Iter;
  CGAL::Polyhedron_3<Gt,I,HDS,A>& ncp = const_cast<CGAL::Polyhedron_3<Gt,I,HDS,A>&>(p);
  return std::make_pair( Iter(ncp.edges_begin()), Iter(ncp.edges_end()) );
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
next_edge_ccw( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor inedge
             , Polyhedron_3<Gt,I,HDS,A> const& 
             )
{
  HalfedgeDS_items_decorator< Polyhedron_3<Gt,I,HDS,A> > D ;
  return D.get_prev(inedge->opposite());
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor
next_edge_cw( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor inedge
            , Polyhedron_3<Gt,I,HDS,A> const& 
            )
{
  return inedge->next()->opposite();
}


//
// Non-const versions
// 

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct halfedge_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> > 
  : CGAL::HDS_halfedge_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >
{
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline std::pair<typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> >::undirected_edge_iterator
                ,typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> >::undirected_edge_iterator 
                >  
undirected_edges( Polyhedron_3<Gt,I,HDS,A>& p )
{
  typedef typename halfedge_graph_traits< Polyhedron_3<Gt,I,HDS,A> >::undirected_edge_iterator Iter;
  CGAL::Polyhedron_3<Gt,I,HDS,A>& ncp = const_cast<CGAL::Polyhedron_3<Gt,I,HDS,A>&>(p);
  return std::make_pair( Iter(ncp.edges_begin()), Iter(ncp.edges_end()) );
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
         , Polyhedron_3<Gt,I,HDS,A>&
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
next_edge_ccw( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor inedge
             , Polyhedron_3<Gt,I,HDS,A>& 
             )
{
  HalfedgeDS_items_decorator< Polyhedron_3<Gt,I,HDS,A> > D ;
  return D.get_prev(inedge->opposite());
}


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor
next_edge_cw( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor inedge
            , Polyhedron_3<Gt,I,HDS,A>& 
            )
{
  return inedge->next()->opposite();
}


} //namespace CGAL

#undef CGAL_HDS_

#endif // CGAL_BOOST_GRAPH_HALFEDGE_GRAPH_TRAITS_POLYHEDRON_3_H
