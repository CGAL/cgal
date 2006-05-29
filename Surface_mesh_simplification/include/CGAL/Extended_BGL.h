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
// $URL: $
// $Id: $
// 
//
// Author(s) : Fernando Caccciola <fernando.cacciola@gmail.com>

#ifndef CGAL_EXTENDED_BGL_H
#define CGAL_EXTENDED_BGL_H

#include <CGAL/Polyhedron_BGL.h>

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
#  define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS
#else
#  define CGAL_HDS_PARAM_ class HDS
#endif

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> > : HDS_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >

CGAL_BEGIN_NAMESPACE
 
template<class HDS>  
typename boost::graph_traits<HDS>::edge_descriptor 
next_edge( typename boost::graph_traits<HDS>::edge_descriptor outedge, HDS& hds)
{
  return outedge->next();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_const_descriptor
next_edge( typename boost::graph_traits<HDS>::edge_const_descriptor outedge, HDS const& hds)
{
  return outedge->next();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_descriptor 
prev_edge( typename boost::graph_traits<HDS>::edge_descriptor outedge, HDS& hds)
{
  CGAL::HalfedgeDS_items_decorator<HDS> D ;
  return D.get_prev(outedge);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_const_descriptor 
prev_edge( typename boost::graph_traits<HDS>::edge_const_descriptor outedge, HDS const& hds)
{
  CGAL::HalfedgeDS_items_decorator<HDS> D ;
  return D.get_prev(outedge);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_descriptor
opposite_edge( typename boost::graph_traits<HDS>::edge_descriptor e, HDS& hds)
{
  return e->opposite();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_const_descriptor
opposite_edge( typename boost::graph_traits<HDS>::edge_const_descriptor e, HDS const& hds)
{
  return e->opposite();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_descriptor 
next_edge_ccw( typename boost::graph_traits<HDS>::edge_descriptor outedge, HDS& hds)
{
  CGAL::HalfedgeDS_items_decorator<HDS> D ;
  return D.get_prev(outedge)->opposite();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_const_descriptor
next_edge_ccw( typename boost::graph_traits<HDS>::edge_const_descriptor outedge, HDS const& hds)
{
  CGAL::HalfedgeDS_items_decorator<HDS> D ;
  return D.get_prev(outedge)->opposite();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_descriptor
next_edge_cw( typename boost::graph_traits<HDS>::edge_descriptor outedge, HDS& hds)
{
  return outedge->opposite()->next();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_const_descriptor
next_edge_cw( typename boost::graph_traits<HDS>::edge_const_descriptor outedge, HDS const& hds)
{
  return outedge->opposite()->next();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_descriptor 
next_edge( typename boost::graph_traits<HDS>::edge_descriptor outedge, HDS& hds)
{
  return outedge->next();
}





template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_descriptor 
out_edge( typename boost::graph_traits<HDS>::vertex_descriptor u, HDS& hds)
{
  HalfedgeDS_items_decorator<HDS> D ;
  return D.get_vertex_edge(u)->opposite();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_const_descriptor
out_edge( typename boost::graph_traits<HDS>::vertex_const_descriptor u, HDS const& hds)
{
  HalfedgeDS_items_decorator<HDS> D ;
  return D.get_vertex_edge(u)->opposite();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_descriptor 
in_edge( typename boost::graph_traits<HDS>::vertex_descriptor u, HDS& hds)
{
  HalfedgeDS_items_decorator<HDS> D ;
  return D.get_vertex_edge(u);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::edge_const_descriptor
in_edge( typename boost::graph_traits<HDS>::vertex_const_descriptor u, HDS const& hds)
{
  HalfedgeDS_items_decorator<HDS> D ;
  return D.get_vertex_edge(u);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::vertex_descriptor
next_vertex( typename boost::graph_traits<HDS>::vertex_descriptor v, HDS& hds)
{
  return boost::target(out_edge(v,hds),hds);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::vertex_const_descriptor
next_vertex( typename boost::graph_traits<HDS>::vertex_const_descriptor v, HDS const& hds)
{
  return boost::target(out_edge(v,hds),hds);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::vertex_descriptor
prev_vertex( typename boost::graph_traits<HDS>::vertex_descriptor v, HDS& hds)
{
  return boost::source(prev_edge(out_edge(v,hds)),hds);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
template<class HDS>  
typename boost::graph_traits<HDS>::vertex_const_descriptor
prev_vertex( typename boost::graph_traits<HDS>::vertex_const_descriptor v, HDS const& hds)
{
  return boost::source(prev_edge(out_edge(v,hds)),hds);
}
        
CGAL_END_NAMESPACE

#undef CGAL_HDS_

#endif // CGAL_POLYHEDRON_EXTENDED_TRAITS_H
