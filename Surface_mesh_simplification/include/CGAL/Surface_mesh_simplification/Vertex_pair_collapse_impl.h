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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_IMPL_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_IMPL_H

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{

template<class M,class D,class C,class V,class S, class I>
VertexPairCollapse<M,D,C,V,S,I>::VertexPairCollapse( TSM&                           aSurface
                                                   , SetCollapseData const&         aSet_collapse_data
                                                   , ParamsToSetCollapseData const* aParamsToSetCollapseData
                                                   , GetCost const&                 aGet_cost
                                                   , GetNewVertexPoint const&       aGet_new_vertex_point
                                                   , ShouldStop const&              aShould_stop
                                                   , VisitorT*                      aVisitor 
                                                   , bool                           aIncludeNonEdgePairs 
                                                   )
  : 
   mSurface (aSurface)
  ,mParamsToSetCollapseData(aParamsToSetCollapseData)
   
  ,Set_collapse_data   (aSet_collapse_data)
  ,Get_cost            (aGet_cost)
  ,Get_new_vertex_point(aGet_new_vertex_point)
  ,Should_stop         (aShould_stop) 
  ,Visitor             (aVisitor)
  
  ,mIncludeNonEdgePairs(aIncludeNonEdgePairs)
{
  CGAL_TSMS_TRACE(0,"VertexPairCollapse of TSM with " << num_undirected_edges(aSurface) << " edges" );
  
  CGAL_TSMS_DEBUG_CODE ( mStep = 0 ; )
  
#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
    vertex_iterator vb, ve ;     
    for ( tie(vb,ve) = vertices(mSurface) ; vb != ve ; ++ vb )  
      CGAL_TSMS_TRACE(2,"V" << (*vb)->ID << xyz_to_string((*vb)->point()) ) ;
      
    undirected_edge_iterator eb, ee ;
    for ( tie(eb,ee) = undirected_edges(mSurface); eb!=ee; ++eb )
      CGAL_TSMS_TRACE(2,"E" << (*eb)->ID << " (V" << source((*eb),mSurface)->ID << "-V" << target((*eb),mSurface)->ID << ")" ) ;
#endif
}

template<class M,class D,class C,class V,class S, class I>
int VertexPairCollapse<M,D,C,V,S,I>::run()
{
  if ( Visitor )
    Visitor->OnStarted(mSurface);
   
  // First collect all candidate edges in a PQ
  Collect(); 
  
  // Then proceed to collapse each edge in turn
  Loop();

  CGAL_TSMS_TRACE(0,"Finished: " << (mInitialPairCount - mCurrentPairCount) << " pairs removed." ) ;

  int r = (int)(mInitialPairCount - mCurrentPairCount) ;
    
  if ( Visitor )
    Visitor->OnFinished(mSurface);
    
  return r ;
}

template<class M,class D,class C,class V,class S, class I>
void VertexPairCollapse<M,D,C,V,S,I>::Collect()
{
  CGAL_TSMS_TRACE(0,"Collecting vertex-pairs...");

  //
  // NOTE: This algorithm requires ALL directed edges to be mapped the corresponding vertex-pair.
  // This must be true whether the pair is in the PQ or not
  //
  
  //
  // Loop over all the edges in the surface putting the corresponding vertex-pairs in the PQ
  //
  
  Equal_3 equal_points = Kernel().equal_3_object();
    
  CGAL_TSMS_DEBUG_CODE( size_type lID = 0 ; ) 
  
  std::size_t lSize = num_undirected_edges(mSurface) ;
  
  mInitialPairCount = mCurrentPairCount = lSize;
  
  mVertexPairArray.reset( new vertex_pair[lSize] ) ;
  
  vertex_pair_ptr lPairPtr = mVertexPairArray.get();
  
  vertex_pair_ptr lPairPtr_End = lPairPtr + lSize ;
  
  undirected_edge_iterator eb, ee ;
  for ( tie(eb,ee) = undirected_edges(mSurface); eb!=ee; ++eb )
  {
    CGAL_assertion( lPairPtr < lPairPtr_End ) ;
    
    edge_descriptor edge = *eb ;
    
    vertex_descriptor s = source(edge,mSurface);
    vertex_descriptor t = target(edge,mSurface);
    
    bool is_s_fixed = is_vertex_fixed(s) ;
    bool is_t_fixed = is_vertex_fixed(t) ;
 
    // All directed edges must be mapped to vertex-pair whether they are put in the PQ or not.
    // Thus, degenerate edges are marked as fixed.
    if ( s == t || equal_points( get_point(s), get_point(t)) )
      is_s_fixed = is_t_fixed = true ;
      
    // The collapse always keeps the vertex 't' and removes vertex 's'.
    // If the 's' is fixed (but not t) then swap the pair to keep (the original) 's' correctly fixed.
    if ( is_s_fixed && !is_t_fixed )
    {
      using std::swap ;
      swap(s,t);
      swap(is_s_fixed,is_t_fixed);
      edge = opposite_edge(edge,mSurface);
    }
    
    CGAL_TSMS_DEBUG_CODE( lPairPtr->id() = lID++ ) ;
    
    Set_collapse_data(lPairPtr->data(),s,t,is_s_fixed,is_t_fixed,edge,mSurface,mParamsToSetCollapseData) ;
    
    lPairPtr->cost() = Get_cost(lPairPtr->data()) ;
    
    set_pair(edge,lPairPtr);
    
    if ( !lPairPtr->is_edge_fixed() )
      insert_in_PQ(lPairPtr);
      
    if ( Visitor )
      Visitor->OnCollected(mSurface,s,t,is_s_fixed,is_t_fixed,edge,lPairPtr->cost(),Get_new_vertex_point(lPairPtr->data()));
    
    CGAL_TSMS_TRACE(2,*lPairPtr);
    
    ++ lPairPtr ;
  }
 
  // TODO: Collect non-edge pairs if requested 
  
  CGAL_TSMS_TRACE(0,"Initial pair count: " << mInitialPairCount ) ;
}

template<class M,class D,class C,class V,class S, class I>
void VertexPairCollapse<M,D,C,V,S,I>::Loop()
{
  CGAL_TSMS_TRACE(0,"Removing pairs...") ;

  //
  // Pops and processes each vertex-pair from the PQ
  //
  vertex_pair_ptr lPair ;
  while ( (lPair = pop_from_PQ()) != 0 )
  {
    CGAL_TSMS_TRACE(3,"Poped " << *lPair ) ;
    
    bool lIsCollapsable = false ;
    
    if ( lPair->cost() != none ) 
    {
      if ( Should_stop(*lPair->cost(),lPair->data(),mInitialPairCount,mCurrentPairCount) )
      {
        if ( Visitor )
          Visitor->OnStopConditionReached(mSurface);
        CGAL_TSMS_TRACE(0,"Stop condition reached with InitialCount=" << mInitialPairCount
                       << " CurrentPairCount=" << mCurrentPairCount
                       << " CurrentPair: " << *lPair
                       );
        break ;
      }

      lIsCollapsable = Is_collapsable(lPair->p(),lPair->q(),lPair->edge()) ;
        
      if ( !lPair->is_edge_fixed() && lIsCollapsable )
      {
        Collapse(lPair);
      }
      else
      {
        CGAL_TSMS_TRACE(1,*lPair << " NOT Collapsable"  );
      }  
    }
    else
    {
      CGAL_TSMS_TRACE(1,*lPair << " uncomputable cost."  );
    }
    
    if ( Visitor )
      Visitor->OnProcessed(mSurface
                          ,lPair->p()
                          ,lPair->q()
                          ,lPair->is_p_fixed()
                          ,lPair->is_q_fixed()
                          ,lPair->edge()
                          ,lPair->cost()
                          ,Get_new_vertex_point(lPair->data())
                          ,lIsCollapsable
                          );
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
bool VertexPairCollapse<M,D,C,V,S,I>::Is_collapsable( vertex_descriptor const& p, vertex_descriptor const& q, edge_descriptor const& p_q )
{
  bool rR = true ;

  CGAL_TSMS_TRACE(3,"Testing collapsabilty of p_q=V" << p->ID << "(%" << p->vertex_degree() << ")"
                 << "->V" << q->ID << "(%" << q->vertex_degree() << ")" 
                 );
                 
  CGAL_TSMS_TRACE(4, "is p_q border:" << is_border(p_q) );
  CGAL_TSMS_TRACE(4, "is q_q border:" << is_border(opposite_edge(p_q,mSurface)) ) ;

  bool is_boundary = is_undirected_edge_a_border(p_q) ;  
  std::size_t min = is_boundary ? 3 : 4 ;
  if ( num_vertices(mSurface) > min )
  {
    out_edge_iterator eb1, ee1 ; 
    out_edge_iterator eb2, ee2 ; 
  
    edge_descriptor q_p = opposite_edge(p_q,mSurface);
    
    vertex_descriptor t = target(next_edge(p_q,mSurface),mSurface);
    vertex_descriptor b = target(next_edge(q_p,mSurface),mSurface);
  
    CGAL_TSMS_TRACE(4,"  t=V" << t->ID << "(%" << t->vertex_degree() << ")" );
    CGAL_TSMS_TRACE(4,"  b=V" << b->ID << "(%" << b->vertex_degree() << ")" );

    // The following loop checks the link condition for p_q.
    // Specifically, that every vertex 'k' adjacent to both 'p and 'q' is a face of the mesh.
    // 
    for ( tie(eb1,ee1) = out_edges(p,mSurface) ; rR && eb1 != ee1 ; ++ eb1 )
    {
      edge_descriptor p_k = *eb1 ;
      
      if ( p_k != p_q )
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
            bool is_face =   ( t == k && !is_border(p_q) )
                          || ( b == k && !is_border(q_p) ) ;
                          
            if ( !is_face )
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
     
  if ( rR && !is_boundary )
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
void VertexPairCollapse<M,D,C,V,S,I>::Collapse( vertex_pair_ptr aPair )
{
  CGAL_TSMS_TRACE(1,"S" << mStep << ". Collapsig " << *aPair ) ;
  
  vertex_descriptor lP = aPair->p();
  vertex_descriptor lQ = aPair->q();
  
  CGAL_assertion( lP != lQ );
  
  // The external function Get_new_vertex_point() is allowed to return an absent point if there is no way to place the vertex
  // satisfying its constrians. In that case the vertex-pair is simply not removed.
  Optional_vertex_point_type lNewVertexPoint = Get_new_vertex_point(aPair->data());
  if ( lNewVertexPoint )
  {
    CGAL_TSMS_TRACE(2,"New vertex point: " << xyz_to_string(*lNewVertexPoint) ) ;

    -- mCurrentPairCount ;
        
    edge_descriptor lEdgePQ = aPair->edge();
    edge_descriptor lEdgeQP = opposite_edge(lEdgePQ,mSurface);
    
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
      vertex_pair_ptr lPairPT = get_pair(lEdgePT) ;
      
      if ( lPairPT->is_in_PQ() )
      {
        CGAL_TSMS_TRACE(2,"Removing VP" << lPairPT->id() << " from PQ" ) ;
        remove_from_PQ(lPairPT) ;
        -- mCurrentPairCount ;
      }
    }
    
    if ( handle_assigned(lEdgeQB) )
    {
      CGAL_TSMS_TRACE(3,"EdgeQB E" << lEdgeQB->ID 
                     << "(V" <<  lEdgeQB->vertex()->ID << "->V" << lEdgeQB->opposite()->vertex()->ID 
                     << ") EdgeBQ E" << lEdgeQB->opposite()->ID 
                     ) ;
      vertex_pair_ptr lPairQB = get_pair(lEdgeQB) ;
      if ( lPairQB->is_in_PQ() )
      {
        CGAL_TSMS_TRACE(2,"Removing VP" << lPairQB->id() << " from PQ") ;
        remove_from_PQ(lPairQB) ;
        -- mCurrentPairCount ;
      }
    }

    if ( Visitor )
      Visitor->OnCollapsed(mSurface,lP,lEdgePQ,lEdgePT,lEdgeQB);
      
    CGAL_TSMS_TRACE(1,"Removing:\n  P-Q: E" << lEdgePQ->ID << "(V" << lP->ID << "->V" << lQ->ID << ")" );
    CGAL_TSMS_TRACE_IF(handle_assigned(lEdgePT),1,"  P-T: E" << lEdgePT->ID << "(V" << lP->ID << "->V" << target(lEdgePT,mSurface)->ID << ")" ) ;
    CGAL_TSMS_TRACE_IF(handle_assigned(lEdgeQB),1,"  Q-B: E" << lEdgeQB->ID << "(V" << lQ->ID << "->V" << target(lEdgeQB,mSurface)->ID << ")" ) ;
    
    // Perform the actuall collapse.
    // This is an external function.
    // It's REQUIRED to remove ONLY 1 vertex (P or Q) and edges PQ,PT and QB (PT and QB are removed if they are not null).
    // All other edges must be kept.
    // All directed edges incident to vertex removed are relink to the vertex kept.
    vertex_descriptor lKeptV = Collapse_triangulation_edge(lEdgePQ,lEdgePT,lEdgeQB,mSurface);
    
    vertex_descriptor lRemovedV = ( lKeptV == lP ? lQ : lP ) ; 
    
    CGAL_TSMS_TRACE(1,"V" << lKeptV->ID << " kept." ) ;
                   
#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
    out_edge_iterator eb1, ee1 ;     
    for ( tie(eb1,ee1) = out_edges(lKeptV,mSurface) ; eb1 != ee1 ; ++ eb1 )  
      CGAL_TSMS_TRACE(2,"  E" << (*eb1)->ID << " V" << lKeptV->ID << "->V" << target(*eb1,mSurface)->ID ) ;
#endif
    
    // Reset the point of placement of the kept vertex.
    vertex_point_t vertex_point ;
    put(vertex_point,mSurface,lKeptV,*lNewVertexPoint) ;

    Update_neighbors(lRemovedV,lKeptV, ( lKeptV == lP ? aPair->is_p_fixed() : aPair->is_q_fixed() ) ) ;
  }
  else
  {
    CGAL_TSMS_TRACE(0,"Unable to calculate new vertex point. Pair " << aPair << " discarded without being removed" ) ;
  }
  
  CGAL_TSMS_DEBUG_CODE ( ++mStep ; )
}

template<class M,class D,class C,class V,class S, class I>
void VertexPairCollapse<M,D,C,V,S,I>::Update_neighbors( vertex_descriptor aRemovedV
                                                      , vertex_descriptor aKeptV
                                                      , bool              aIsKeptVFixed
                                                      )
{
  CGAL_TSMS_TRACE(3,"Updating cost of neighboring edges..." ) ;

  //
  // (A) Collect all pairs to update its cost: all those around each vertex adjacent to the vertex kept
  //
  
  typedef std::set<vertex_pair_ptr> vertex_pair_set ;
  
  vertex_pair_set lToUpdate ;
  
  // (A.1) Loop around all vertices adjacent to the vertex kept
  in_edge_iterator eb1, ee1 ; 
  for ( tie(eb1,ee1) = in_edges(aKeptV,mSurface) ; eb1 != ee1 ; ++ eb1 )
  {
    edge_descriptor edge1 = *eb1 ;
    
    vertex_descriptor adj_k = source(edge1,mSurface);
    
    // (A.2) Loop around all edges incident on each adjacent vertex
    in_edge_iterator eb2, ee2 ; 
    for ( tie(eb2,ee2) = in_edges(adj_k,mSurface) ; eb2 != ee2 ; ++ eb2 )
    {
      edge_descriptor edge2 = *eb2 ;
      
      vertex_pair_ptr lPair = get_pair(edge2);
      
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
      
      CGAL_TSMS_TRACE(4,"Inedge around V" << adj_k->ID << " E" << edge2->ID << " Opposite E" << edge2->opposite()->ID 
                     << " V" <<  edge2->opposite()->vertex()->ID << "->V" << edge2->vertex()->ID 
                     << "\nPair:" << *lPair
                     ) ;
                   
    
      // Only those pairs still in the PQ are update.
      // The mark is used because in the way we loop here the same pair is found many times.
      if ( lPair->is_in_PQ() && lToUpdate.find(lPair) == lToUpdate.end() )
      {
        CGAL_TSMS_TRACE(4,"Pair registered for updating.") ;
        lToUpdate.insert(lPair) ;
      }  
    } 
  }  
  
  //
  // (B) Proceed to update the costs.
  //
  
  for ( typename vertex_pair_set::iterator it = lToUpdate.begin(), eit = lToUpdate.end() ; it != eit ; ++ it )
  {
    vertex_pair_ptr lPair = *it ;

    Set_collapse_data(lPair->data()
                     ,lPair->p()
                     ,lPair->q()
                     ,lPair->is_p_fixed()     
                     ,lPair->is_q_fixed()
                     ,lPair->edge()
                     ,lPair->surface()
                     ,mParamsToSetCollapseData
                     ) ;
    lPair->cost() = Get_cost(lPair->data());                 
    
    CGAL_TSMS_TRACE(3,"Updated cost of " << *lPair) ;
    
    update_in_PQ(lPair);
  }
    
}


} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_IMPL_H //
// EOF //
 
