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
// $URL$
// $Id$
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_OPERATOR_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_OPERATOR_H 1

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{


//
// collapse_triangulation_edge euler operator (topological operation only).
//
// This operator collapses the edge p-q replacing it with one single vertex connected to the link of p-q.
//
// The net effect of the operator is equivalent to removing one of the vertices and re-triangulating the resulting hole.
//
// The actual collapse operation only removes up to 2 faces, 3 edges and 1 vertex, as follows.
//
// The top face, p-q-t, to the left of p->q, if any, is removed.
// The bottom face, q-p-b, to the left of q->p (right of p->q), if any, is removed.
//
// If there is a top-left face, that is, if the edge p->t is NOT a border edge, the top face is removed by joining it
// with the top-left face (this creates a 4 sided face: p-q-t-l).
// If p->t IS a border edge then the top face is simply erased.
//
// If there is a bottom-right face, that is, if the edge q->b is NOT a border edge, the bottom face is removed by joining it
// with the bottom-right face (this creates a 4 sided face: q-p-b-r)
// If q->b IS a border edge then the bottom face is simply erased.
//
// If there is a top face edge p->t is removed. 
// If there is a bottom face edge q->b is removed. 
// 
// One of the vertices (p or q) is removed. 
// If there is no top-left face so the top face is directly erased AND there is no bottom face,
// that face erasure automatically removes vertex p. Likewise, if there is no bottom-right face 
// so the bottom face is directly erased AND there is no top face, that erasure automatically 
// removes vertex q.
// Directly erasing the top/bottom faces when the opposite face exists does not removes any vertex
// automatically, in which case vertex p is removed by joining p->q
//
// NOTES:
//
// This operator can only be called by collapsable edges (which satisfy the link condition), hence,
// there must exist at least the top face or a bottom face, and, if there is no top-left face there is a top-right face
// and likewise, if there is no bottom-right face there is a bottom-left face. (IOW vertices t/b must have degree >=3
// unless the top/bottom faces do not exists)
//
// The operator doesn't join the top face with the top-left face and the bottom-face with the bottom-left face
// (or both to the right) because if the top-left and bottom-left (or both right) faces are themselve adjacent,
// the first joint would introduce a degree 2 vertex.
// That is why the operator alternates left and right to join the top and bottom faces.
//
// PARAMETERS:
//  pq : the edge to collapse
//  pt : the edge shared between the top face and the top-left face (if any). If there is no top face this parameter is a null handle.
//  qb : the edge shared between the bottom face and the bottom-right face. If there is no bottom face this parameter is a null handle.
//
// RETURN VALUE: A handle to the vertex that IS NOT removed.
//
template<class DS_>
struct Collapse_triangulation_edge
{
  typedef DS_ DS ;
  
  typedef typename boost::graph_traits<DS>::edge_descriptor   edge_descriptor ;
  typedef typename boost::graph_traits<DS>::vertex_descriptor vertex_descriptor ;
  
  vertex_descriptor operator() ( edge_descriptor const& pq
                               , edge_descriptor const& pt
                               , edge_descriptor const& qb
                               , DS&                    aDS 
                               ) const
  {
    edge_descriptor null ;
    
    bool lTopFaceExists         = pt != null ;
    bool lBottomFaceExists      = qb != null ;
    bool lTopLeftFaceExists     = lTopFaceExists    && !pt->is_border() ;
    bool lBottomRightFaceExists = lBottomFaceExists && !qb->is_border() ;
    
    CGAL_precondition( !lTopFaceExists    || (lTopFaceExists    && ( pt->vertex()->vertex_degree() > 2 ) ) ) ;
    CGAL_precondition( !lBottomFaceExists || (lBottomFaceExists && ( qb->vertex()->vertex_degree() > 2 ) ) ) ;
    
    vertex_descriptor q = pq->vertex();
    vertex_descriptor p = pq->opposite()->vertex();
    
    CGAL_TSMS_TRACE(3, "Collapsing p-q E" << pq->ID << " (V" << p->ID << "->V" << q->ID << ")" ) ;
    
    bool lP_Erased = false, lQ_Erased = false ;
    
    if ( lTopFaceExists )
    { 
      CGAL_precondition( !pt->opposite()->is_border() ) ; // p-q-t is a face of the mesh
      if ( lTopLeftFaceExists )
      {
        CGAL_TSMS_TRACE(3, "Removing p-t E" << pt->ID << " (V" << p->ID << "->V" << pt->vertex()->ID << ") by joining top-left face" ) ;
        
        aDS.join_facet (pt);
      }
      else
      {
        CGAL_TSMS_TRACE(3, "Removing p-t E" << pt->ID << " (V" << p->ID << "->V" << pt->vertex()->ID << ") by erasing top face" ) ;
        
        aDS.erase_facet(pt->opposite());
        
        if ( !lBottomFaceExists )
        {
          CGAL_TSMS_TRACE(3, "Bottom face doesn't exist so vertex P already removed" ) ;
          lP_Erased = true ;
        }  
      } 
    }
    
    if ( lBottomFaceExists )
    {   
      CGAL_precondition( !qb->opposite()->is_border() ) ; // p-q-b is a face of the mesh
      if ( lBottomRightFaceExists )
      {
        CGAL_TSMS_TRACE(3, "Removing q-b E" << qb->ID << " (V" << q->ID << "->V" << qb->vertex()->ID << ") by joining bottom-right face" ) ;
        aDS.join_facet (qb);
      }
      else
      {
        CGAL_TSMS_TRACE(3, "Removing q-b E" << qb->ID << " (V" << q->ID << "->V" << qb->vertex()->ID << ") by erasing bottom face" ) ;
        
        aDS.erase_facet(qb->opposite());
        
        if ( !lTopFaceExists )
        {
          CGAL_TSMS_TRACE(3, "Top face doesn't exist so vertex Q already removed" ) ;
          lQ_Erased = true ;
        }  
      }
    }

    CGAL_assertion( !lP_Erased || !lQ_Erased ) ;
    
    if ( !lP_Erased && !lQ_Erased )
    {
      CGAL_TSMS_TRACE(3, "Removing vertex P by joining pQ" ) ;
      aDS.join_vertex(pq);
      lP_Erased = true ;
    }    
    
    return lP_Erased ? q : p ;
  }
  
} ;

} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_OPERATOR_H //
// EOF //
 
