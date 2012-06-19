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

#ifndef CGAL_BOOST_GRAPH_GRAPH_TRAITS_POLYHEDRON_3_H
#define CGAL_BOOST_GRAPH_GRAPH_TRAITS_POLYHEDRON_3_H

#include <CGAL/boost/graph/graph_traits_HalfedgeDS.h>

#include <CGAL/Polyhedron_3.h>

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

namespace boost
{ 

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >
   : CGAL::HDS_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >
{};


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const >
   : CGAL::HDS_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> > // See NOTE above!
{};

 
template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertices_size_type
num_vertices(const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  return p.size_of_vertices();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edges_size_type 
num_edges(const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  return p.size_of_halfedges() ;
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::degree_size_type
degree(typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor v, const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  return v->vertex_degree() * 2 ;
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::degree_size_type
out_degree(typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor v, const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  return v->vertex_degree();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::degree_size_type
in_degree(typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor v,  const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  return v->vertex_degree();
}


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor
source(typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor e, const CGAL::Polyhedron_3<Gt,I,HDS,A> & )
{
  return e->opposite()->vertex();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor
target(typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor e, const CGAL::Polyhedron_3<Gt,I,HDS,A> & )
{
  return e->vertex();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline std::pair<typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_iterator
                ,typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_iterator 
                >  
vertices( const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  typedef typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_iterator Iter;
  CGAL::Polyhedron_3<Gt,I,HDS,A>& ncp = const_cast<CGAL::Polyhedron_3<Gt,I,HDS,A>&>(p);
  return std::make_pair( Iter(ncp.vertices_begin()), Iter(ncp.vertices_end()) );
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline std::pair<typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edge_iterator
                ,typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edge_iterator 
                >  
edges( const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  typedef typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edge_iterator Iter;
  CGAL::Polyhedron_3<Gt,I,HDS,A>& ncp = const_cast<CGAL::Polyhedron_3<Gt,I,HDS,A>&>(p);
  return std::make_pair( Iter(ncp.halfedges_begin()), Iter(ncp.halfedges_end()) );
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline std::pair<typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::in_edge_iterator
                ,typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::in_edge_iterator
                >  
in_edges( typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor u, const CGAL::Polyhedron_3<Gt,I,HDS,A>& g)
{
  typename CGAL::Polyhedron_3<Gt,I,HDS,A>::Halfedge_around_vertex_circulator ec = u->vertex_begin();
  typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edges_size_type in_deg = in_degree(u,g);
  typedef typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::in_edge_iterator Iter;
  return std::make_pair( Iter(ec), Iter(ec,in_deg) );
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline std::pair<typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::out_edge_iterator
                ,typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::out_edge_iterator 
                >  
out_edges( typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor u, const CGAL::Polyhedron_3<Gt,I,HDS,A>& g)
{
  typename CGAL::Polyhedron_3<Gt,I,HDS,A>::Halfedge_around_vertex_circulator ec = u->vertex_begin();
  typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edges_size_type out_deg = out_degree(u,g);
  typedef typename graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::out_edge_iterator Iter;
  return std::make_pair( Iter(ec), Iter(ec,out_deg) );
}
        
} // namespace boost

#undef CGAL_HDS_PARAM_

#endif // CGAL_BOOST_GRAPH_GRAPH_TRAITS_POLYHEDRON_3_H
