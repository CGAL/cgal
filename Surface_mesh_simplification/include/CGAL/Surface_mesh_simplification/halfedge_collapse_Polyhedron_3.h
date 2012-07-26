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

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor
halfedge_collapse( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor const& pq
                 , Polyhedron_3<Gt,I,HDS,A>& aSurface
                 )     
{
  typedef Polyhedron_3<Gt,I,HDS,A> Surface ;
  
  typedef typename boost::graph_traits<Surface>::vertex_descriptor vertex_descriptor ;
  typedef typename boost::graph_traits<Surface>::edge_descriptor   edge_descriptor ;
  
  edge_descriptor qp = opposite_edge(pq,aSurface);
  edge_descriptor pt = opposite_edge(prev_edge(pq,aSurface),aSurface);
  edge_descriptor qb = opposite_edge(prev_edge(qp,aSurface),aSurface);
  
  bool lTopFaceExists         = !pq->is_border() ;
  bool lBottomFaceExists      = !qp->is_border() ;
  bool lTopLeftFaceExists     = lTopFaceExists    && !pt->is_border() ;
  bool lBottomRightFaceExists = lBottomFaceExists && !qb->is_border() ;
  
  CGAL_precondition( !lTopFaceExists    || (lTopFaceExists    && ( pt->vertex()->vertex_degree() > 2 ) ) ) ;
  CGAL_precondition( !lBottomFaceExists || (lBottomFaceExists && ( qb->vertex()->vertex_degree() > 2 ) ) ) ;
  
  vertex_descriptor q = pq->vertex();
  vertex_descriptor p = pq->opposite()->vertex();
  
  CGAL_ECMS_TRACE(3, "Collapsing p-q E" << pq->id() << " (V" << p->id() << "->V" << q->id() << ")" ) ;
  
  bool lP_Erased = false, lQ_Erased = false ;
  
  if ( lTopFaceExists )
  { 
    CGAL_precondition( !pt->opposite()->is_border() ) ; // p-q-t is a face of the mesh
    if ( lTopLeftFaceExists )
    {
      CGAL_ECMS_TRACE(3, "Removing p-t E" << pt->id() << " (V" << p->id() << "->V" << pt->vertex()->id() << ") by joining top-left face" ) ;
      
      aSurface.join_facet (pt);
    }
    else
    {
      CGAL_ECMS_TRACE(3, "Removing p-t E" << pt->id() << " (V" << p->id() << "->V" << pt->vertex()->id() << ") by erasing top face" ) ;
      
      aSurface.erase_facet(pt->opposite());
      
      if ( !lBottomFaceExists )
      {
        CGAL_ECMS_TRACE(3, "Bottom face doesn't exist so vertex P already removed" ) ;
        lP_Erased = true ;
      }  
    } 
  }
  
  if ( lBottomFaceExists )
  {   
    CGAL_precondition( !qb->opposite()->is_border() ) ; // p-q-b is a face of the mesh
    if ( lBottomRightFaceExists )
    {
      CGAL_ECMS_TRACE(3, "Removing q-b E" << qb->id() << " (V" << q->id() << "->V" << qb->vertex()->id() << ") by joining bottom-right face" ) ;
      aSurface.join_facet (qb);
    }
    else
    {
      CGAL_ECMS_TRACE(3, "Removing q-b E" << qb->id() << " (V" << q->id() << "->V" << qb->vertex()->id() << ") by erasing bottom face" ) ;
      
      aSurface.erase_facet(qb->opposite());
      
      if ( !lTopFaceExists )
      {
        CGAL_ECMS_TRACE(3, "Top face doesn't exist so vertex Q already removed" ) ;
        lQ_Erased = true ;
      }  
    }
  }

  CGAL_assertion( !lP_Erased || !lQ_Erased ) ;
  
  if ( !lP_Erased && !lQ_Erased )
  {
    CGAL_ECMS_TRACE(3, "Removing vertex P by joining pQ" ) ;
    aSurface.join_vertex(pq);
    lP_Erased = true ;
  }    
  
  return lP_Erased ? q : p ;
} 

} // namespace Surface_mesh_simplification

} //namespace CGAL

#undef CGAL_HDS_PARAM

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_TRIANGULATION_EDGE_POLYHEDRON_3_H
// EOF //
 
