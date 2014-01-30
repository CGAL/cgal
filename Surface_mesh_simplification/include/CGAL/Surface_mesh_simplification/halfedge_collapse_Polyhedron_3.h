// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_TRIANGULATION_EDGE_POLYHEDRON_3_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_TRIANGULATION_EDGE_POLYHEDRON_3_H

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Polyhedron_3.h>

#  define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS

namespace CGAL {

namespace Surface_mesh_simplification
{
/*
Function responsible for contracting an edge while respecting constrained edges

Notations used in the following function:
Top=TopFace
Btm=BottomFace

          t
        /   \
       /     \
      /  Top  \
     p -------- q
      \  Btm  /
       \     /
        \   /
          b

Prerequisites:
If Top exists, amongst p-t and t-q only one is constrained
If Btm exists, amongst p-b and q-b only one is constrained
p-q is not constrained
*/
template<class Gt, class I, CGAL_HDS_PARAM_, class A, class EdgeIsConstrainedMap>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor
halfedge_collapse( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor const& pq
                 , Polyhedron_3<Gt,I,HDS,A>& aSurface
                 , EdgeIsConstrainedMap Edge_is_constrained_map
                 )
{
  CGAL_assertion( !get(Edge_is_constrained_map,pq) );
  typedef Polyhedron_3<Gt,I,HDS,A> Surface ;

  typedef typename boost::graph_traits<Surface>::vertex_descriptor vertex_descriptor ;
  typedef typename boost::graph_traits<Surface>::edge_descriptor   edge_descriptor ;

  edge_descriptor qp = opposite_edge(pq,aSurface);
  edge_descriptor pt = opposite_edge(prev_edge(pq,aSurface),aSurface);
  edge_descriptor qb = opposite_edge(prev_edge(qp,aSurface),aSurface);
  edge_descriptor tq = pq->next()->opposite();
  edge_descriptor bp = qp->next()->opposite();

  bool lTopFaceExists         = !pq->is_border() ;
  bool lBottomFaceExists      = !qp->is_border() ;

  CGAL_precondition( !lTopFaceExists    || (lTopFaceExists    && ( pt->vertex()->vertex_degree() > 2 ) ) ) ;
  CGAL_precondition( !lBottomFaceExists || (lBottomFaceExists && ( qb->vertex()->vertex_degree() > 2 ) ) ) ;

  vertex_descriptor q = pq->vertex();
  vertex_descriptor p = pq->opposite()->vertex();

  CGAL_ECMS_TRACE(3, "Collapsing p-q E" << pq->id() << " (V" << p->id() << "->V" << q->id() << ")" ) ;

  //used to collect edges to remove from the surface
  edge_descriptor edges_to_erase[2];
  edge_descriptor* edges_to_erase_ptr=edges_to_erase;

  // If the top facet exists, we need to choose one out of the two edges which one disappears:
  //   p-t if it is not constrained and t-q otherwise
  if ( lTopFaceExists )
  {
    CGAL_precondition( !pt->opposite()->is_border() ) ; // p-q-t is a face of the mesh
    if ( !get(Edge_is_constrained_map,pt) )
    {
      CGAL_ECMS_TRACE(3, "Removing p-t E" << pt->id() << " (V" << p->id() << "->V" << pt->vertex()->id()) ;
      *edges_to_erase_ptr++=pt;
    }
    else
    {
      CGAL_ECMS_TRACE(3, "Removing t-q E" << pt->id() << " (V" << pt->vertex()->id() << "->V" << q->id() ) ;
      CGAL_assertion( !get(Edge_is_constrained_map,tq) );
      *edges_to_erase_ptr++=tq;
    }
  }

  // If the bottom facet exists, we need to choose one out of the two edges which one disappears:
  //   q-b if it is not constrained and b-p otherwise
  if ( lBottomFaceExists )
  {
    if ( !get(Edge_is_constrained_map,qb) )
    {
      CGAL_ECMS_TRACE(3, "Removing q-b E" << qb->id() << " (V" << q->id() << "->V" << qb->vertex()->id() ) ;
      *edges_to_erase_ptr++=qb;
    }
    else{
      CGAL_ECMS_TRACE(3, "Removing b-p E" << qb->id() << " (V" << qb->vertex()->id() << "->V" << p->id() ) ;
      CGAL_assertion( !get(Edge_is_constrained_map,bp) );
      *edges_to_erase_ptr++=bp;
    }
  }

  if (lTopFaceExists && lBottomFaceExists)
  {
    if ( edges_to_erase[0]->facet()==edges_to_erase[1]->facet()
         && !edges_to_erase[0]->is_border() )
    {
      // the vertex is of valence 3 and we simply need to remove the vertex
      // and its indicent edges
      bool lP_Erased=false;
      edge_descriptor edge =
        edges_to_erase[0]->next()==edges_to_erase[1]?
          edges_to_erase[0]:edges_to_erase[1];
      if (edge->vertex()==p)
        lP_Erased=true;
      aSurface.erase_center_vertex(edge);
      return lP_Erased? q : p;
    }
    else
    {
      if (!edges_to_erase[0]->is_border())
        aSurface.join_facet(edges_to_erase[0]);
      else
        aSurface.erase_facet(edges_to_erase[0]->opposite());
      if (!edges_to_erase[1]->is_border())
        aSurface.join_facet(edges_to_erase[1]);
      else
        aSurface.erase_facet(edges_to_erase[1]->opposite());
      aSurface.join_vertex(pq);
      return q;
    }
  }
  else
  {
      if (lTopFaceExists)
      {
        if (!edges_to_erase[0]->is_border()){
          aSurface.join_facet(edges_to_erase[0]);
          aSurface.join_vertex(pq);
          return q;
        }
        bool lQ_Erased=pq->next()->opposite()->is_border();
        aSurface.erase_facet(edges_to_erase[0]->opposite());
        return lQ_Erased?p:q;
      }

      CGAL_assertion(lBottomFaceExists);
      if (!edges_to_erase[0]->is_border()){
        aSurface.join_facet(edges_to_erase[0]);
        aSurface.join_vertex(qp);
        return p;
      }
      bool lP_Erased=qp->next()->opposite()->is_border();
      aSurface.erase_facet(edges_to_erase[0]->opposite());
      return lP_Erased?q:p;
  };
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor
halfedge_collapse( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor const& pq
                 , Polyhedron_3<Gt,I,HDS,A>& aSurface
                 )
{
  return halfedge_collapse( pq,
                            aSurface,
                            No_constrained_edge_map<Polyhedron_3<Gt,I,HDS,A> >() );
}

} // namespace Surface_mesh_simplification

} //namespace CGAL

#undef CGAL_HDS_PARAM

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_TRIANGULATION_EDGE_POLYHEDRON_3_H
// EOF //

