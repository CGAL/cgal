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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_COLLAPSE_IMPL_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_COLLAPSE_IMPL_H

CGAL_BEGIN_NAMESPACE

namespace Surface_mesh_simplification 
{

template<class M,class SP,class VPM, class VFM, class EIM,class EBM, class SC, class CF,class PF,class CP, class PP,class V>
EdgeCollapse<M,SP,VPM,VFM,EIM,EBM,SC,CF,PF,CP,PP,V>::EdgeCollapse( ECM&                    aSurface
                                                                 , ShouldStop       const& aShould_stop
                                                                 , VertexPointMap   const& aVertex_point_map 
                                                                 , VertexIsFixedMap const& aVertex_is_fixed_map 
                                                                 , EdgeIndexMap     const& aEdge_index_map 
                                                                 , EdgeIsBorderMap  const& aEdge_is_border_map 
                                                                 , SetCache         const& aSet_cache
                                                                 , GetCost          const& aGet_cost
                                                                 , GetPlacement     const& aGet_placement
                                                                 , CostParams       const* aCostParams
                                                                 , PlacementParams  const* aPlacementParams
                                                                 , VisitorT*               aVisitor 
                                                                 )
  : 
   mSurface           (aSurface)
  ,Should_stop        (aShould_stop) 
  ,Vertex_point_map   (aVertex_point_map)
  ,Vertex_is_fixed_map(aVertex_is_fixed_map)
  ,Edge_index_map     (aEdge_index_map)
  ,Edge_is_border_map (aEdge_is_border_map)
  ,Set_cache          (aSet_cache)
  ,Get_cost           (aGet_cost)
  ,Get_placement      (aGet_placement)
  ,mCostParams        (aCostParams)
  ,mPlacementParams   (aPlacementParams)
  ,Visitor            (aVisitor)
  
{
  CGAL_ECMS_TRACE(0,"EdgeCollapse of ECM with " << (num_edges(aSurface)/2) << " edges" ); 
  
  CGAL_ECMS_DEBUG_CODE ( mStep = 0 ; )
  
#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
    vertex_iterator vb, ve ;     
    for ( tie(vb,ve) = vertices(mSurface) ; vb != ve ; ++ vb )  
      CGAL_ECMS_TRACE(2, vertex_to_string(*vb) ) ;
      
    undirected_edge_iterator eb, ee ;
    for ( tie(eb,ee) = undirected_edges(mSurface); eb!=ee; ++eb )
      CGAL_ECMS_TRACE(2, edge_to_string(*eb) ) ;
#endif
}

template<class M,class SP,class VPM, class VFM, class EIM,class EBM, class SC, class CF,class PF,class CP, class PP,class V>
int EdgeCollapse<M,SP,VPM,VFM,EIM,EBM,SC,CF,PF,CP,PP,V>::run()
{
  if ( Visitor )
    Visitor->OnStarted(mSurface);
   
  // First collect all candidate edges in a PQ
  Collect(); 
  
  // Then proceed to collapse each edge in turn
  Loop();

  CGAL_ECMS_TRACE(0,"Finished: " << (mInitialEdgeCount - mCurrentEdgeCount) << " edges removed." ) ;

  int r = (int)(mInitialEdgeCount - mCurrentEdgeCount) ;
    
  if ( Visitor )
    Visitor->OnFinished(mSurface);
    
  return r ;
}

template<class M,class SP,class VPM, class VFM, class EIM,class EBM, class SC, class CF,class PF,class CP, class PP,class V>
void EdgeCollapse<M,SP,VPM,VFM,EIM,EBM,SC,CF,PF,CP,PP,V>::Collect()
{
  CGAL_ECMS_TRACE(0,"Collecting edges...");

  //
  // Loop over all the _undirected_ edges in the surface putting them in the PQ
  //
  
  Equal_3 equal_points = Kernel().equal_3_object();
    
  size_type lSize = num_edges(mSurface) / 2 ;
  
  mInitialEdgeCount = mCurrentEdgeCount = lSize;
  
  mEdgeDataArray.reset( new Edge_data[lSize] ) ;
  
  mPQ.reset( new PQ (lSize, Compare_cost(this), Undirected_edge_id(this) ) ) ;
  
  std::size_t id = 0 ;
  
  undirected_edge_iterator eb, ee ;
  for ( tie(eb,ee) = undirected_edges(mSurface); eb!=ee; ++eb )
  {
    edge_descriptor lEdge = *eb ;
  
    CGAL_assertion( get_directed_edge_id(lEdge) == id ) ;
    CGAL_assertion( get_directed_edge_id(opposite_edge(lEdge,mSurface)) == id+1 ) ;
      
    vertex_descriptor p,q ;
    tie(p,q) = get_vertices(lEdge);
    
    bool lIsFixed = is_vertex_fixed(p) || is_vertex_fixed(q) ;
 
    if ( p == q || equal_points( get_point(p), get_point(q)) )
      lIsFixed = true ;
  
    // But in the case of fixed edges the edge data is left default constructed
    if ( !lIsFixed )
    {
      Edge_data& lData = get_data(lEdge);
      Set_cache(lData.cache(),lEdge,mSurface,mCostParams,mPlacementParams) ;
      insert_in_PQ(lEdge,lData);
    }
      
    if ( Visitor )
      Visitor->OnCollected(lEdge,lIsFixed,mSurface);
    
    CGAL_ECMS_TRACE(2,edge_to_string(lEdge));
    
    id += 2 ;
  }
 
  CGAL_ECMS_TRACE(0,"Initial edge count: " << mInitialEdgeCount ) ;
}

template<class M,class SP,class VPM, class VFM, class EIM,class EBM, class SC, class CF,class PF,class CP, class PP,class V>
void EdgeCollapse<M,SP,VPM,VFM,EIM,EBM,SC,CF,PF,CP,PP,V>::Loop()
{
  CGAL_ECMS_TRACE(0,"Collapsing edges...") ;

  //
  // Pops and processes each edge from the PQ
  //
  optional<edge_descriptor> lEdge ;
  while ( (lEdge = pop_from_PQ()) )
  {
    CGAL_ECMS_TRACE(3,"Poped " << edge_to_string(*lEdge) ) ;
    
    Optional_cost_type lCost = get_cost(*lEdge);
    
    if ( Visitor )
      Visitor->OnSelected(*lEdge,mSurface,lCost,mInitialEdgeCount,mCurrentEdgeCount);
      
    if ( lCost != none ) 
    {
      if ( Should_stop(*lCost,*lEdge,mInitialEdgeCount,mCurrentEdgeCount) )
      {
        if ( Visitor )
          Visitor->OnStopConditionReached(mSurface);
          
        CGAL_ECMS_TRACE(0,"Stop condition reached with InitialEdgeCount=" << mInitialEdgeCount
                       << " CurrentEdgeCount=" << mCurrentEdgeCount
                       << " Current Edge: " << edge_to_string(*lEdge)
                       );
        break ;
      }
        
      if ( Is_collapsable(*lEdge) )
      {
        Collapse(*lEdge);
      }
      else
      {
        if ( Visitor )
          Visitor->OnNonCollapsable(*lEdge,mSurface);
          
        CGAL_ECMS_TRACE(1,edge_to_string(*lEdge) << " NOT Collapsable"  );
      }  
    }
    else
    {
      CGAL_ECMS_TRACE(1,edge_to_string(*lEdge) << " uncomputable cost."  );
    }
    
  }
}

template<class M,class SP,class VPM, class VFM, class EIM,class EBM, class SC, class CF,class PF,class CP, class PP,class V>
bool EdgeCollapse<M,SP,VPM,VFM,EIM,EBM,SC,CF,PF,CP,PP,V>::is_border( const_vertex_descriptor const& aV ) const
{
  bool rR = false ;
  
  const_in_edge_iterator eb, ee ; 
  for ( tie(eb,ee) = in_edges(aV,mSurface) ; eb != ee ; ++ eb )
  {
    const_edge_descriptor lEdge = *eb ;
    if ( is_undirected_edge_a_border(lEdge) )
    {
      rR = true ;
      break ;
    }
  }  
    
  return rR ;  
}

// Some edges are NOT collapsable: doing so would break the topological consistency of the mesh.
// This function returns true if a edge 'p->q' can be collapsed.
//
// An edge p->q can be collapsed iff it satisfies the "link condition"
// (as described in the "Mesh Optimization" article of Hoppe et al (1993))
//
// The link conidition is as follows: for every vertex 'k' adjacent to both 'p and 'q', "p,k,q" is a facet of the mesh.
//
template<class M,class SP,class VPM, class VFM, class EIM,class EBM, class SC, class CF,class PF,class CP, class PP,class V>
bool EdgeCollapse<M,SP,VPM,VFM,EIM,EBM,SC,CF,PF,CP,PP,V>::Is_collapsable( edge_descriptor const& aEdgePQ )
{
  bool rR = true ;

  vertex_descriptor p,q ; tie(p,q) = get_vertices(aEdgePQ);
  
  CGAL_ECMS_TRACE(3,"Testing collapsabilty of p_q=V" << p->ID << "(%" << p->vertex_degree() << ")"
                 << "->V" << q->ID << "(%" << q->vertex_degree() << ")" 
                 );
                 
  CGAL_ECMS_TRACE(4, "is p_q border:" << is_border(aEdgePQ) );
  CGAL_ECMS_TRACE(4, "is q_q border:" << is_border(opposite_edge(aEdgePQ,mSurface)) ) ;

  bool lIsBoundary = is_undirected_edge_a_border(aEdgePQ) ;  
  std::size_t min = lIsBoundary ? 3 : 4 ;
  if ( num_vertices(mSurface) > min )
  {
    out_edge_iterator eb1, ee1 ; 
    out_edge_iterator eb2, ee2 ; 
  
    edge_descriptor lEdgeQP = opposite_edge(aEdgePQ,mSurface);
    
    vertex_descriptor t = target(next_edge(aEdgePQ,mSurface),mSurface);
    vertex_descriptor b = target(next_edge(lEdgeQP,mSurface),mSurface);
  
    CGAL_ECMS_TRACE(4,"  t=V" << t->ID << "(%" << t->vertex_degree() << ")" );
    CGAL_ECMS_TRACE(4,"  b=V" << b->ID << "(%" << b->vertex_degree() << ")" );

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
              CGAL_ECMS_TRACE(3,"  k=V" << k->ID << " IS NOT in a face with p-q. NON-COLLAPSABLE edge." ) ;
              rR = false ;
              break ;
            }  
            else 
            {
              CGAL_ECMS_TRACE(4,"  k=V" << k->ID << " is in a face with p-q") ;
            }
          }
        }  
      }
    }   
  }
  else
  {
    rR = false ;
    CGAL_ECMS_TRACE(3,"  Surface is irreducible. NON-COLLAPSABLE edge." ) ;
  }
     
  if ( rR && !lIsBoundary )
  {
    if ( is_border(p) && is_border(q) )
    {
      rR = false ;
      CGAL_ECMS_TRACE(3,"  both p and q are boundary vertices but p-q is not. NON-COLLAPSABLE edge." ) ;
    }  
  }
  
  return rR ;
}

template<class M,class SP,class VPM, class VFM, class EIM,class EBM, class SC, class CF,class PF,class CP, class PP,class V>
void EdgeCollapse<M,SP,VPM,VFM,EIM,EBM,SC,CF,PF,CP,PP,V>::Collapse( edge_descriptor const& aEdgePQ )
{
  CGAL_ECMS_TRACE(1,"S" << mStep << ". Collapsig " << edge_to_string(aEdgePQ) ) ;
  
  vertex_descriptor lP, lQ ; tie(lP,lQ) = get_vertices(aEdgePQ);
  
  CGAL_assertion( lP != lQ );

  vertex_descriptor rResult ;
    
  // The external function Get_new_vertex_point() is allowed to return an absent point if there is no way to place the vertex
  // satisfying its constrians. In that case the vertex-pair is simply not removed.
  Optional_placement_type lPlacement = get_placement(aEdgePQ);
  
  if ( Visitor )
    Visitor->OnCollapsing(aEdgePQ,mSurface,lPlacement);
    
  if ( lPlacement )
  {
    CGAL_ECMS_TRACE(2,"New vertex point: " << xyz_to_string(*lPlacement) ) ;

    -- mCurrentEdgeCount ;
        
    edge_descriptor lEdgeQP = opposite_edge(aEdgePQ,mSurface);
    
    // The collapse of P-Q removes the top and bottom facets, if any.
    // Edges P-T and P-Q, which are in the top and bottom facets (if they exist), are used by the collapse operator to remove them.
    
    // Edges P-T and P-Q are defined only if the top/bottom facets exists
    
    edge_descriptor lEdgePT, lEdgeQB ;
    
    if ( !is_border(aEdgePQ) ) // Exists top facet
      lEdgePT = opposite_edge(prev_edge(aEdgePQ,mSurface),mSurface);
    
    if ( !is_border(lEdgeQP) ) // Exists bottom facet
      lEdgeQB = opposite_edge(prev_edge(lEdgeQP,mSurface),mSurface);
    
    CGAL_ECMS_TRACE(3,"EdgePQ E" << aEdgePQ->ID 
                   << "(V" <<  aEdgePQ->vertex()->ID << "->V" << aEdgePQ->opposite()->vertex()->ID 
                   << ") EdgeQP E" << aEdgePQ->opposite()->ID 
                   ) ;
                   

    // If the top/bottom facets exists, they are removed and the edges P-T and Q-B along with them.
    // In that case their corresponding pairs must be pop off the queue
    
    if ( handle_assigned(lEdgePT) )
    {
      CGAL_ECMS_TRACE(3,"EdgePT E" << lEdgePT->ID 
                     << "(V" <<  lEdgePT->vertex()->ID << "->V" << lEdgePT->opposite()->vertex()->ID 
                     << ") EdgeTP E" << lEdgePT->opposite()->ID 
                     ) ;
                     
      Edge_data& lDataPT = get_data(lEdgePT) ;
      if ( lDataPT.is_in_PQ() )
      {
        CGAL_ECMS_TRACE(2,"Removing E" << lEdgePT->ID << " from PQ" ) ;
        remove_from_PQ(lEdgePT,lDataPT) ;
        -- mCurrentEdgeCount ;
      }
    }
    
    if ( handle_assigned(lEdgeQB) )
    {
      CGAL_ECMS_TRACE(3,"EdgeQB E" << lEdgeQB->ID 
                     << "(V" <<  lEdgeQB->vertex()->ID << "->V" << lEdgeQB->opposite()->vertex()->ID 
                     << ") EdgeBQ E" << lEdgeQB->opposite()->ID 
                     ) ;
                     
      Edge_data& lDataQB = get_data(lEdgeQB) ;
      if ( lDataQB.is_in_PQ() )
      {
        CGAL_ECMS_TRACE(2,"Removing E" << lEdgeQB->ID << " from PQ") ;
        remove_from_PQ(lEdgeQB,lDataQB) ;
        -- mCurrentEdgeCount ;
      }
    }

    CGAL_ECMS_TRACE(1,"Removing:\n  P-Q: E" << aEdgePQ->ID << "(V" << lP->ID << "->V" << lQ->ID << ")" );
    CGAL_ECMS_TRACE_IF(handle_assigned(lEdgePT),1,"  P-T: E" << lEdgePT->ID << "(V" << lP->ID << "->V" << target(lEdgePT,mSurface)->ID << ")" ) ;
    CGAL_ECMS_TRACE_IF(handle_assigned(lEdgeQB),1,"  Q-B: E" << lEdgeQB->ID << "(V" << lQ->ID << "->V" << target(lEdgeQB,mSurface)->ID << ")" ) ;
    
      
    // Perform the actuall collapse.
    // This is an external function.
    // It's REQUIRED to remove ONLY 1 vertex (P or Q) and edges PQ,PT and QB (PT and QB are removed if they are not null).
    // All other edges must be kept.
    // All directed edges incident to vertex removed are relink to the vertex kept.
    rResult = collapse_triangulation_edge(aEdgePQ,mSurface);
    
    CGAL_ECMS_TRACE(1,"V" << rResult->ID << " kept." ) ;
                   
#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
    out_edge_iterator eb1, ee1 ;     
    for ( tie(eb1,ee1) = out_edges(rResult,mSurface) ; eb1 != ee1 ; ++ eb1 )  
      CGAL_ECMS_TRACE(2, edge_to_string(*eb1) ) ;
#endif
    
    put(vertex_point,mSurface,rResult,*lPlacement) ;

    Update_neighbors(rResult) ;
  }
  else
  {
    CGAL_ECMS_TRACE(0,"Unable to calculate new vertex point. E" << aEdgePQ->ID << " discarded without being collapsed" ) ;
  }
  
  CGAL_ECMS_DEBUG_CODE ( ++mStep ; )
}

template<class M,class SP,class VPM, class VFM, class EIM,class EBM, class SC, class CF,class PF,class CP, class PP,class V>
void EdgeCollapse<M,SP,VPM,VFM,EIM,EBM,SC,CF,PF,CP,PP,V>::Update_neighbors( vertex_descriptor const& aKeptV )
{
  CGAL_ECMS_TRACE(3,"Updating cost of neighboring edges..." ) ;

  //
  // (A) Collect all edges to update its cost: all those around each vertex adjacent to the vertex kept
  //  
  
  typedef std::set<edge_descriptor,Compare_id> edges ;
  
  edges lToUpdate(Compare_id(this)) ;
  
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
      
      Edge_data& lData2 = get_data(lEdge2);
      CGAL_ECMS_TRACE(4,"Inedge around V" << lAdj_k->ID << edge_to_string(lEdge2) ) ;
    
      // Only those edges still in the PQ _and_ not already collected are updated.
      if ( lData2.is_in_PQ() && lToUpdate.find(lEdge2) == lToUpdate.end() )
        lToUpdate.insert(lEdge2) ;
    } 
  }  
  
  //
  // (B) Proceed to update the costs.
  //
  
  for ( typename edges::iterator it = lToUpdate.begin(), eit = lToUpdate.end() ; it != eit ; ++ it )
  {
    edge_descriptor lEdge = *it;
    
    Edge_data& lData = get_data(lEdge);
    
    Set_cache(lData.cache(),lEdge,mSurface,mCostParams,mPlacementParams) ;
    
    CGAL_ECMS_TRACE(3, edge_to_string(lEdge) << " updated in the PQ") ;
    
    update_in_PQ(lEdge,lData);
  }
    
}


} // namespace Surface_mesh_simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_COLLAPSE_IMPL_H //
// EOF //
 
