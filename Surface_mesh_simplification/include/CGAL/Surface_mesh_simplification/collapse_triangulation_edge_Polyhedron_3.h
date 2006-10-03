 // Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Surface_mesh_simplification/Collapse_operator.h $
// $Id: Collapse_operator.h 33792 2006-08-25 21:51:01Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_TRIANGULATION_EDGE_POLYHEDRON_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_TRIANGULATION_EDGE_POLYHEDRON_H

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Polyhedron_3.h>

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
#  define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS
#else
#  define CGAL_HDS_PARAM_ class HDS
#endif

CGAL_BEGIN_NAMESPACE

namespace Surface_mesh_simplification
{

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor
collapse_triangulation_edge( typename boost::graph_traits< Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor const& pq
                           , Polyhedron_3<Gt,I,HDS,A>& aSurface
                           )     
{
  edge_descriptor qp = opposite_edge(pq,aSurface);
  
  bool lTopFaceExists         = !pq->is_border() ;
  bool lBottomFaceExists      = !qp->is_border() ;
  bool lTopLeftFaceExists     = lTopFaceExists    && !pt->is_border() ;
  bool lBottomRightFaceExists = lBottomFaceExists && !qb->is_border() ;
  
  edge_descriptor pq = opposite_edge(prev_edge(pq,aSurface),aSurface);
  edge_descriptor qb = opposite_edge(prev_edge(qp,aSurface),aSurface);
  
  CGAL_precondition( !lTopFaceExists    || (lTopFaceExists    && ( pt->vertex()->vertex_degree() > 2 ) ) ) ;
  CGAL_precondition( !lBottomFaceExists || (lBottomFaceExists && ( qb->vertex()->vertex_degree() > 2 ) ) ) ;
  
  vertex_descriptor q = pq->vertex();
  vertex_descriptor p = pq->opposite()->vertex();
  
  CGAL_ECMS_TRACE(3, "Collapsing p-q E" << pq->ID << " (V" << p->ID << "->V" << q->ID << ")" ) ;
  
  bool lP_Erased = false, lQ_Erased = false ;
  
  if ( lTopFaceExists )
  { 
    CGAL_precondition( !pt->opposite()->is_border() ) ; // p-q-t is a face of the mesh
    if ( lTopLeftFaceExists )
    {
      CGAL_ECMS_TRACE(3, "Removing p-t E" << pt->ID << " (V" << p->ID << "->V" << pt->vertex()->ID << ") by joining top-left face" ) ;
      
      aSurface.join_facet (pt);
    }
    else
    {
      CGAL_ECMS_TRACE(3, "Removing p-t E" << pt->ID << " (V" << p->ID << "->V" << pt->vertex()->ID << ") by erasing top face" ) ;
      
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
      CGAL_ECMS_TRACE(3, "Removing q-b E" << qb->ID << " (V" << q->ID << "->V" << qb->vertex()->ID << ") by joining bottom-right face" ) ;
      aSurface.join_facet (qb);
    }
    else
    {
      CGAL_ECMS_TRACE(3, "Removing q-b E" << qb->ID << " (V" << q->ID << "->V" << qb->vertex()->ID << ") by erasing bottom face" ) ;
      
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

} // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#undef CGAL_HDS_

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_TRIANGULATION_EDGE_POLYHEDRON_H
// EOF //
 
