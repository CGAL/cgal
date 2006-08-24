// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_I_EDGE_COLLAPSE_IMPL_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_I_EDGE_COLLAPSE_IMPL_H

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{

template<class M,class D,class C,class V,class S, class I>
EdgeCollapse<M,D,C,V,S,I>::EdgeCollapse( TSM&                           aSurface
                                       , SetCollapseData const&         aSet_collapse_data
                                       , ParamsToSetCollapseData const* aParamsToSetCollapseData
                                       , GetCost const&                 aGet_cost
                                       , GetNewVertexPoint const&       aGet_new_vertex_point
                                       , ShouldStop const&              aShould_stop
                                       , VisitorT*                      aVisitor 
                                       )
  : 
   mSurface (aSurface)
  ,mParamsToSetCollapseData(aParamsToSetCollapseData)
   
  ,Set_collapse_data   (aSet_collapse_data)
  ,Get_cost            (aGet_cost)
  ,Get_new_vertex_point(aGet_new_vertex_point)
  ,Should_stop         (aShould_stop) 
  ,Visitor             (aVisitor)
  
  ,mPQ( Compare_cost(this) )
  
{
  CGAL_TSMS_TRACE(0,"EdgeCollapse of TSM with " << num_undirected_edges(aSurface) << " edges" );
  
  CGAL_TSMS_DEBUG_CODE ( mStep = 0 ; )
  
#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
    vertex_iterator vb, ve ;     
    for ( tie(vb,ve) = vertices(mSurface) ; vb != ve ; ++ vb )  
      CGAL_TSMS_TRACE(2, vertex_to_string(*vb) ) ;
      
    undirected_edge_iterator eb, ee ;
    for ( tie(eb,ee) = undirected_edges(mSurface); eb!=ee; ++eb )
      CGAL_TSMS_TRACE(2, edge_to_string(*eb) ) ;
#endif
}

template<class M,class D,class C,class V,class S, class I>
int EdgeCollapse<M,D,C,V,S,I>::run()
{
  if ( Visitor )
    Visitor->OnStarted(mSurface);
   
  // First collect all candidate edges in a PQ
  Collect(); 
  
  // Then proceed to collapse each edge in turn
  Loop();

  CGAL_TSMS_TRACE(0,"Finished: " << (mInitialEdgeCount - mCurrentEdgeCount) << " edges removed." ) ;

  int r = (int)(mInitialEdgeCount - mCurrentEdgeCount) ;
    
  if ( Visitor )
    Visitor->OnFinished(mSurface);
    
  return r ;
}

template<class M,class D,class C,class V,class S, class I>
void EdgeCollapse<M,D,C,V,S,I>::Collect()
{
  CGAL_TSMS_TRACE(0,"Collecting edges...");

  //
  // Loop over all the edges in the surface putting them in the PQ
  //
  
  Equal_3 equal_points = Kernel().equal_3_object();
    
  CGAL_TSMS_DEBUG_CODE( size_type lID = 0 ; ) 
  
  std::size_t lSize = num_undirected_edges(mSurface) ;
  
  mInitialEdgeCount = mCurrentEdgeCount = lSize;
  
  mEdgeDataArray.reset( new Edge_data[lSize] ) ;
  
  Edge_data_ptr lData = mEdgeDataArray.get();
  
  Edge_data_ptr lData_End = lData + lSize ;
  
  undirected_edge_iterator eb, ee ;
  for ( tie(eb,ee) = undirected_edges(mSurface); eb!=ee; ++eb )
  {
    CGAL_assertion( lData < lData_End ) ;
    
    edge_descriptor lEdge = *eb ;
    
    vertex_descriptor p,q ;
    tie(p,q) = get_vertices(lEdge);
    
    bool lIsFixed = is_vertex_fixed(p) || is_vertex_fixed(q) ;
 
    if ( p == q || equal_points( get_point(p), get_point(q)) )
      is_fixed = true ;
      
    // For simplicity, ALL edges, fixed or not, are associated with an edge data.
    set_data(lEdge,lData);
    
    // But in the case of fixed edges the edge data is left default constructed
    if ( !lIsFixed )
    {
      Set_collapse_data(lData->data(),lEdge,mSurface,mParamsToSetCollapseData) ;
      insert_in_PQ(lEdge,lData);
    }
      
    if ( Visitor )
      Visitor->OnCollected(lEdge,lIsFixed,mSurface);
    
    CGAL_TSMS_TRACE(2,edge_to_string(lEdge));
    
    ++ lData ;
  }
 
  CGAL_TSMS_TRACE(0,"Initial edge count: " << mInitialEdgeCount ) ;
}

template<class M,class D,class C,class V,class S, class I>
void EdgeCollapse<M,D,C,V,S,I>::Loop()
{
  CGAL_TSMS_TRACE(0,"Collapsing edges...") ;

  //
  // Pops and processes each edge from the PQ
  //
  edge_descriptor lEdge ;
  while ( handle_assigned(lEdge = pop_from_PQ()) )
  {
    CGAL_TSMS_TRACE(3,"Poped " << edge_to_string(lEdge) ) ;
    
    bool lIsCollapsable = false ;
    
    vertex_descritor lVertex ;
    
    Optional_cost_type lCost = get_cost(lEdge);
    
    if ( lCost != none ) 
    {
      if ( Should_stop(*lCost,lEdge,mInitialEdgeCount,mCurrentEdgeCount) )
      {
        if ( Visitor )
          Visitor->OnStopConditionReached(mSurface);
          
        CGAL_TSMS_TRACE(0,"Stop condition reached with InitialEdgeCount=" << mInitialEdgeCount
                       << " CurrentEdgeCount=" << mCurrentEdgeCount
                       << " Current Edge: " << edge_to_string(lEdge)
                       );
        break ;
      }
        
      if ( Is_collapsable(lEdge) )
      {
        lVertex= Collapse(lEdge);
      }
      else
      {
        CGAL_TSMS_TRACE(1,edge_to_string(lEdge) << " NOT Collapsable"  );
      }  
    }
    else
    {
      CGAL_TSMS_TRACE(1,edge_to_string(lEdge) << " uncomputable cost."  );
    }
    
    if ( Visitor )
      Visitor->OnProcessed(lEdge,mSurface,lCost,lVertex);
  }
}

// Some edges are NOT collapsable: doing so would break the topological consistency of the mesh.
// This function returns true if a edge 'p->q' can be collapsed.
//
// An edge p->q can be collapsed iff it satisfies the "link condition"
// (as described in the "Mesh Optimization" article of Hoppe et al (1993))
//
// The link conidition is as follows: for every vertex 'k' adjacent to both 'p and 'q', "p,k,q" is a facet of the mesh.
//
template<class M,class D,class C,class V,class S, class I>
bool EdgeCollapse<M,D,C,V,S,I>::Is_collapsable( edge_descriptor const& aEdgePQ )
{
  bool rR = true ;

  vertex_descriptor p,q ; tie(p,q) = get_vertices(aEdgePQ);
  
  CGAL_TSMS_TRACE(3,"Testing collapsabilty of p_q=V" << p->ID << "(%" << p->vertex_degree() << ")"
                 << "->V" << q->ID << "(%" << q->vertex_degree() << ")" 
                 );
                 
  CGAL_TSMS_TRACE(4, "is p_q border:" << is_border(aEdgePQ) );
  CGAL_TSMS_TRACE(4, "is q_q border:" << is_border(opposite_edge(aEdgePQ,mSurface)) ) ;

  bool lIsBoundary = is_undirected_edge_a_border(aEdgePQ) ;  
  std::size_t min = lIsBoundary ? 3 : 4 ;
  if ( num_vertices(mSurface) > min )
  {
    out_edge_iterator eb1, ee1 ; 
    out_edge_iterator eb2, ee2 ; 
  
    edge_descriptor lEdgeQP = opposite_edge(aEdgePQ,mSurface);
    
    vertex_descriptor t = target(next_edge(aEdgePQ,mSurface),mSurface);
    vertex_descriptor b = target(next_edge(lEdgeQP,mSurface),mSurface);
  
    CGAL_TSMS_TRACE(4,"  t=V" << t->ID << "(%" << t->vertex_degree() << ")" );
    CGAL_TSMS_TRACE(4,"  b=V" << b->ID << "(%" << b->vertex_degree() << ")" );

    // The following loop checks the link condition for aEdgePQ.
    // Specifically, that every vertex 'k' adjacent to both 'p and 'q' is a face of the mesh.
    // 
    for ( tie(eb1,ee1) = out_edges(p,mSurface) ; rR && eb1 != ee1 ; ++ eb1 )
    {
      edge_descriptor p_k = *eb1 ;
      
      if ( p_k != aEdgePQ )
      {
        vertex_descriptor k = target(p_k,mSurface);
        
        for ( tie(eb2,ee2) = out_edges(k,mSurface) ; rR && eb2 != ee2 ; ++ eb2 )
        {
          edge_descriptor k_l = *eb2 ;
  
          if ( target(k_l,mSurface) == q )
          {
            // At this point we know p-q-k are connected and we need to determine if this triangle is a face of the mesh.
            //
            // Since the mesh is known to be triangular there are at most two faces sharing the edge p-q.
            //
            // If p->q is NOT a border edge, the top face is p->q->t where t is target(next(p->q))
            // If q->p is NOT a border edge, the bottom face is q->p->b where b is target(next(q->p))
            //
            // If k is either t or b then p-q-k *might* be a face of the mesh. It won't be if k==t but p->q is border
            // or k==b but q->b is a border (because in that case even though there exists triangles p->q->t (or q->p->b)
            // they are holes, not faces)
            // 
            bool lIsFace =   ( t == k && !is_border(aEdgePQ) )
                          || ( b == k && !is_border(lEdgeQP) ) ;
                          
            if ( !lIsFace )
            {
              CGAL_TSMS_TRACE(3,"  k=V" << k->ID << " IS NOT in a face with p-q. NON-COLLAPSABLE edge." ) ;
              rR = false ;
            }  
            else 
            {
              CGAL_TSMS_TRACE(4,"  k=V" << k->ID << " is in a face with p-q") ;
            }
          }
        }  
      }
    }   
  }
  else
  {
    rR = false ;
    CGAL_TSMS_TRACE(3,"  Surface is irreducible. NON-COLLAPSABLE edge." ) ;
  }
     
  if ( rR && !lIsBoundary )
  {
    if ( is_border(p) && is_border(q) )
    {
      rR = false ;
      CGAL_TSMS_TRACE(3,"  both p and q are boundary vertices but p-q is not. NON-COLLAPSABLE edge." ) ;
    }  
  }
  
  return rR ;
}

template<class M,class D,class C,class V,class S, class I>
EdgeCollapse<M,D,C,V,S,I>::vertex_descriptor EdgeCollapse<M,D,C,V,S,I>::Collapse( edge_descriptor const& aEdgePQ )
{
  CGAL_TSMS_TRACE(1,"S" << mStep << ". Collapsig " << edge_to_string(aEdgePQ) ) ;
  
  vertex_descriptor lP, lQ ; tie(lP,lQ) = get_vertices(aEdgePQ);
  
  CGAL_assertion( lP != lQ );

  vertex_descriptor rResult ;
    
  // The external function Get_new_vertex_point() is allowed to return an absent point if there is no way to place the vertex
  // satisfying its constrians. In that case the vertex-pair is simply not removed.
  Optional_vertex_point_type lNewVertexPoint = get_new_vertex_point(aEdgePQ);
  if ( lNewVertexPoint )
  {
    CGAL_TSMS_TRACE(2,"New vertex point: " << xyz_to_string(*lNewVertexPoint) ) ;

    -- mCurrentEdgeCount ;
        
    edge_descriptor lEdgeQP = opposite_edge(aEdgePQ,mSurface);
    
    // The collapse of P-Q removes the top and bottom facets, if any.
    // Edges P-T and P-Q, which are in the top and bottom facets (if they exist), are used by the collapse operator to remove them.
    
    // Edges P-T and P-Q are defined only if the top/bottom facets exists
    
    edge_descriptor lEdgePT, lEdgeQB ;
    
    if ( !is_border(lEdgePQ) ) // Exists top facet
      lEdgePT = opposite_edge(prev_edge(lEdgePQ,mSurface),mSurface);
    
    if ( !is_border(lEdgeQP) ) // Exists bottom facet
      lEdgeQB = opposite_edge(prev_edge(lEdgeQP,mSurface),mSurface);
    
    CGAL_TSMS_TRACE(3,"EdgePQ E" << lEdgePQ->ID 
                   << "(V" <<  lEdgePQ->vertex()->ID << "->V" << lEdgePQ->opposite()->vertex()->ID 
                   << ") EdgeQP E" << lEdgePQ->opposite()->ID 
                   ) ;
                   

    // If the top/bottom facets exists, they are removed and the edges P-T and Q-B along with them.
    // In that case their corresponding pairs must be pop off the queue
    
    if ( handle_assigned(lEdgePT) )
    {
      CGAL_TSMS_TRACE(3,"EdgePT E" << lEdgePT->ID 
                     << "(V" <<  lEdgePT->vertex()->ID << "->V" << lEdgePT->opposite()->vertex()->ID 
                     << ") EdgeTP E" << lEdgePT->opposite()->ID 
                     ) ;
      Edge_data_ptr lDataPT = get_data(lEdgePT) ;
      
      if ( lDataPT->is_in_PQ() )
      {
        CGAL_TSMS_TRACE(2,"Removing E" << lEdgePT->ID << " from PQ" ) ;
        remove_from_PQ(lDataPT) ;
        -- mCurrentEdgeCount ;
      }
    }
    
    if ( handle_assigned(lEdgeQB) )
    {
      CGAL_TSMS_TRACE(3,"EdgeQB E" << lEdgeQB->ID 
                     << "(V" <<  lEdgeQB->vertex()->ID << "->V" << lEdgeQB->opposite()->vertex()->ID 
                     << ") EdgeBQ E" << lEdgeQB->opposite()->ID 
                     ) ;
      Edge_data_ptr lDataQB = get_pair(lDataQB) ;
      if ( lDataQB->is_in_PQ() )
      {
        CGAL_TSMS_TRACE(2,"Removing E" << lEdgeQB->ID << " from PQ") ;
        remove_from_PQ(lDataQB) ;
        -- mCurrentEdgeCount ;
      }
    }

    CGAL_TSMS_TRACE(1,"Removing:\n  P-Q: E" << lEdgePQ->ID << "(V" << lP->ID << "->V" << lQ->ID << ")" );
    CGAL_TSMS_TRACE_IF(handle_assigned(lEdgePT),1,"  P-T: E" << lEdgePT->ID << "(V" << lP->ID << "->V" << target(lEdgePT,mSurface)->ID << ")" ) ;
    CGAL_TSMS_TRACE_IF(handle_assigned(lEdgeQB),1,"  Q-B: E" << lEdgeQB->ID << "(V" << lQ->ID << "->V" << target(lEdgeQB,mSurface)->ID << ")" ) ;
    
    // Perform the actuall collapse.
    // This is an external function.
    // It's REQUIRED to remove ONLY 1 vertex (P or Q) and edges PQ,PT and QB (PT and QB are removed if they are not null).
    // All other edges must be kept.
    // All directed edges incident to vertex removed are relink to the vertex kept.
    rResult = Collapse_triangulation_edge(lEdgePQ,lEdgePT,lEdgeQB,mSurface);
    
    CGAL_TSMS_TRACE(1,"V" << rResult->ID << " kept." ) ;
                   
#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
    out_edge_iterator eb1, ee1 ;     
    for ( tie(eb1,ee1) = out_edges(rResult,mSurface) ; eb1 != ee1 ; ++ eb1 )  
      CGAL_TSMS_TRACE(2, edge_to_string(*eb1) ) ;
#endif
    
    // Reset the point of placement of the kept vertex.
    vertex_point_t vertex_point ;
    put(vertex_point,mSurface,rResult,*lNewVertexPoint) ;

    Update_neighbors(rResult) ;
  }
  else
  {
    CGAL_TSMS_TRACE(0,"Unable to calculate new vertex point. E" << aEdgePQ->ID << " discarded without being collapsed" ) ;
  }
  
  CGAL_TSMS_DEBUG_CODE ( ++mStep ; )
  
  return rResult ;
}

template<class M,class D,class C,class V,class S, class I>
void EdgeCollapse<M,D,C,V,S,I>::Update_neighbors( vertex_descriptor const& aKeptV )
{
  CGAL_TSMS_TRACE(3,"Updating cost of neighboring edges..." ) ;

  //
  // (A) Collect all edges to update its cost: all those around each vertex adjacent to the vertex kept
  //
  
  typedef std::vector<edge_descriptor> edges ;
  
  edges lToUpdate ;
  
  // (A.1) Loop around all vertices adjacent to the vertex kept
  in_edge_iterator eb1, ee1 ; 
  for ( tie(eb1,ee1) = in_edges(aKeptV,mSurface) ; eb1 != ee1 ; ++ eb1 )
  {
    edge_descriptor lEdge1 = *eb1 ;
    
    vertex_descriptor lAdj_k = source(lEdge1,mSurface);
    
    // (A.2) Loop around all edges incident on each adjacent vertex
    in_edge_iterator eb2, ee2 ; 
    for ( tie(eb2,ee2) = in_edges(lAdj_k,mSurface) ; eb2 != ee2 ; ++ eb2 )
    {
      edge_descriptor lEdge2 = *eb2 ;
      
      Edge_data_ptr lData2 = get_data(lEdge2);

      /*
      vertex_descriptor p,q ;       
      Edge_data_ptr lPair = get_pair(edge2);
      
      if ( lPair->p() == aRemovedV )
      {
        CGAL_TSMS_TRACE(4,"Replacing lPair->p() with V" << aKeptV->id() << "(Q)" ) ;
        
        CGAL_assertion( aKeptV != lPair->q() ) ;
        
        lPair->data().set(aKeptV
                         ,lPair->q()
                         ,aIsKeptVFixed 
                         ,lPair->is_q_fixed()
                         ,source(edge2,mSurface) == aKeptV ? edge2 : opposite_edge(edge2,mSurface)
                         ,mSurface
                         );
      }
      else if ( lPair->q() == aRemovedV )
      {
        CGAL_TSMS_TRACE(4,"Replacing lPair->q() with V" << aKeptV->id() << " (Q)" ) ;
        
        CGAL_assertion( aKeptV != lPair->p() ) ;
        
        lPair->data().set(lPair->p()
                         ,aKeptV
                         ,lPair->is_p_fixed()
                         ,aIsKeptVFixed 
                         ,target(edge2,mSurface) == aKeptV ? edge2 : opposite_edge(edge2,mSurface)
                         ,mSurface
                         );
      }
      */
      
      CGAL_TSMS_TRACE(4,"Inedge around V" << lAdj_k->ID << edge_to_string(lEdge2) ) ;
    
      // Only those edges still in the PQ are updated.
      // The mark is used because in the way we loop here the same edge is found many times.
      if ( lData2->is_in_PQ() && lToUpdate.find(lEdge2) == lToUpdate.end() )
        lToUpdate.push_back(lPair) ;
    } 
  }  
  
  //
  // (B) Proceed to update the costs.
  //
  
  for ( typename edges::iterator it = lToUpdate.begin(), eit = lToUpdate.end() ; it != eit ; ++ it )
  {
    edge_descriptor lEdge = *it;
    
    Edge_data_ptr lData = get_data(lData);

    Set_collapse_data(lData->data(),lEdge,mSurface,mParamsToSetCollapseData) ;
    
    CGAL_TSMS_TRACE(3, edge_to_string(lEdge) << " updated in the PQ") ;
    
    update_in_PQ(lData);
  }
    
}


} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_I_EDGE_COLLAPSE_IMPL_H //
// EOF //
 
