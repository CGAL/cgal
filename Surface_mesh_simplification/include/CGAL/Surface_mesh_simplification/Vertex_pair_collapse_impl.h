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
  
}

template<class M,class D,class C,class V,class S, class I>
int VertexPairCollapse<M,D,C,V,S,I>::run()
{
  if ( Visitor )
    Visitor->OnStarted(mSurface);
    
  if ( num_vertices(mSurface) > 4 )
  {
    // First collect all candidate edges in a PQ
    Collect(); 
    
    // Then proceed to collapse each edge in turn
    Loop();
  
    CGAL_TSMS_TRACE(0,"Finished: " << (mInitialPairCount - mCurrentPairCount) << " pairs removed." ) ;
  }
  else
  {
    if ( Visitor )
      Visitor->OnThetrahedronReached(mSurface);
    
    mInitialPairCount = mCurrentPairCount = num_undirected_edges(mSurface) ;
    
    CGAL_TSMS_TRACE(0,"A thetrahedron cannot be simplified any further.");
  }

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
    
    CGAL_TSMS_TRACE(3,*lPairPtr);
    
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
      if ( num_vertices(mSurface) <= 4 )
      {
        if ( Visitor )
          Visitor->OnThetrahedronReached(mSurface);
        CGAL_TSMS_TRACE(0,"Thetrahedron reached.");
        break ;
      }
      
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
template<class M,class D,class C,class V,class S, class I>
bool VertexPairCollapse<M,D,C,V,S,I>::Is_collapsable( vertex_descriptor const& p, vertex_descriptor const& q, edge_descriptor const& p_q )
{
  bool rR = true ;
  
  edge_descriptor q_p = opposite_edge(p_q,mSurface);
  
  edge_descriptor p_t = next_edge_ccw(p_q,mSurface);      
  edge_descriptor p_b = next_edge_cw (p_q,mSurface);      
  edge_descriptor q_t = next_edge_cw (q_p,mSurface);      
  edge_descriptor q_b = next_edge_ccw(q_p,mSurface);      

  CGAL_TSMS_TRACE(5, "Testing link condition:"
                  << "\np_q: V" << source(p_q,mSurface)->id() << "->V" << target(p_q,mSurface)->id()
                  << "\nq_p: V" << source(q_p,mSurface)->id() << "->V" << target(q_p,mSurface)->id()
                  << "\np_t: V" << source(p_t,mSurface)->id() << "->V" << target(p_t,mSurface)->id()
                  << "\np_b: V" << source(p_b,mSurface)->id() << "->V" << target(p_b,mSurface)->id()
                  << "\nq_t: V" << source(q_t,mSurface)->id() << "->V" << target(q_t,mSurface)->id()
                  << "\nq_n: V" << source(q_b,mSurface)->id() << "->V" << target(q_b,mSurface)->id()
                  );
                    
  
  // degree(p) and degree(q) >= 3
  if (    target(p_t,mSurface) != q
       && target(p_b,mSurface) != q
       && target(q_t,mSurface) != p
       && target(q_b,mSurface) != p
     )
  {
    // link('p') .intersection. link('q') == link('p_q') (that is, exactly {'t','b'})
    if (    target(p_t,mSurface) == target(q_t,mSurface) 
         && target(p_b,mSurface) == target(q_b,mSurface) 
       )
    {
      edge_descriptor p_tn = p_t ;
      edge_descriptor p_bn = p_b ;
      edge_descriptor q_tn = q_t ;
      edge_descriptor q_bn = q_b ;
      
      do
      {
        p_tn = next_edge_ccw(p_tn,mSurface);      
        p_bn = next_edge_cw (p_bn,mSurface);      
        q_tn = next_edge_cw (q_tn,mSurface);      
        q_bn = next_edge_ccw(q_bn,mSurface);      
        
        CGAL_TSMS_TRACE(5, "\n  p_tn: V" << source(p_tn,mSurface)->id() << "->V" << target(p_tn,mSurface)->id()
                        << "\n  q_tn: V" << source(q_tn,mSurface)->id() << "->V" << target(q_tn,mSurface)->id()
                        << "\n  p_bn: V" << source(p_bn,mSurface)->id() << "->V" << target(p_bn,mSurface)->id()
                        << "\n  q_bn: V" << source(q_bn,mSurface)->id() << "->V" << target(q_bn,mSurface)->id()
                       );
                       
        if (    target(p_tn,mSurface) == target(q_tn,mSurface)
             || target(p_bn,mSurface) == target(q_bn,mSurface)
           )
        {
          rR = false ;
          break ;
        }     
      }
      while ( p_tn != p_bn && q_tn != q_bn ) ;
    }   
    else
    {
      CGAL_TSMS_TRACE(5, "t and b vertices not shared by p and q");
      rR = false ;
    }
  }   
  else
  {
    CGAL_TSMS_TRACE(5, "degree(p) or degree(q) < 3");
    rR = false ;
  }
  
  return rR ;
}

template<class M,class D,class C,class V,class S, class I>
void VertexPairCollapse<M,D,C,V,S,I>::Collapse( vertex_pair_ptr aPair )
{
  CGAL_TSMS_TRACE(1,"Collapsig " << *aPair ) ;
  
  vertex_descriptor lP = aPair->p();
  vertex_descriptor lQ = aPair->q();
  
  CGAL_assertion( lP != lQ );
  
  // The external function Get_new_vertex_point() is allowed to return an absent point if there is no way to place the vertex
  // satisfying its constrians. In that case the vertex-pair is simply not removed.
  Optional_vertex_point_type lNewVertexPoint = Get_new_vertex_point(aPair->data());
  if ( lNewVertexPoint )
  {
    CGAL_TSMS_TRACE(2,"New vertex point: " << xyz_to_string(*lNewVertexPoint) ) ;
        
    // The actual collapse of edge PQ merges the top and bottom facets with its left and right adjacents resp, then
    // joins P and Q.
    
    edge_descriptor lEdgePQ = aPair->edge();
    
    edge_descriptor lEdgeQP = opposite_edge(lEdgePQ,mSurface);
    
    edge_descriptor lEdgePT = next_edge_ccw(lEdgePQ,mSurface);
    
    edge_descriptor lEdgeQB = next_edge_ccw(lEdgeQP,mSurface);
    
    CGAL_TSMS_TRACE(3,"EdgePQ E" << lEdgePQ->ID << " Opposite EdgePQ E" << lEdgePQ->opposite()->ID 
                   << " V" <<  lEdgePQ->opposite()->vertex()->ID << "->V" << lEdgePQ->vertex()->ID 
                   ) ;
    CGAL_TSMS_TRACE(3,"EdgePT E" << lEdgePT->ID << " Opposite EdgePT E" << lEdgePT->opposite()->ID 
                   << " V" <<  lEdgePT->opposite()->vertex()->ID << "->V" << lEdgePT->vertex()->ID 
                   ) ;

    // The collapse will also remove QB so it must be pop off the PQ as well
    vertex_pair_ptr lPairQB = get_pair(lEdgeQB) ;
    vertex_pair_ptr lPairPT = get_pair(lEdgePT) ;
    
    if ( lPairPT->is_in_PQ() )
    {
      CGAL_TSMS_TRACE(2,"Removing VP" << lPairPT->id() << " from PQ" ) ;
      remove_from_PQ(lPairPT) ;
    }
    if ( lPairQB->is_in_PQ() )
    {
      CGAL_TSMS_TRACE(2,"Removing VP" << lPairQB->id() << " from PQ") ;
      remove_from_PQ(lPairQB) ;
    }
    
    CGAL_TSMS_TRACE(2,"Removing from surface V" << lP->ID 
                   << " E" << lEdgePQ->ID
                   << " E" << lEdgePT->ID
                   << " E" << lEdgeQB->ID
                   );

    if ( Visitor )
      Visitor->OnCollapsed(mSurface,lP,lEdgePQ,lEdgePT,lEdgeQB);
      
    //
    // Perform the actuall collapse.
    // This is an external function.
    // It's REQUIRED to remove ONLY vertex P and edges PQ,PT and QB, keeping all other edges and
    // to relink all directed edges incident to P with Q.
    // The algorithm is based on the stability of the remaining edges.
    Collapse_triangulation_edge(lEdgePQ,lEdgePT,lEdgeQB,mSurface);
    
    // Reset the point of placement of Q (the vertex that "replaces" the collapsed edge)
    vertex_point_t vertex_point ;
    put(vertex_point,mSurface,lQ,*lNewVertexPoint) ;

    Update_neighbors(aPair);
            
    mCurrentPairCount -= 3 ;
  }
  else
  {
    CGAL_TSMS_TRACE(0,"Unable to calculate new vertex point. Pair " << aPair << " discarded without being removed" ) ;
  }
}

template<class M,class D,class C,class V,class S, class I>
void VertexPairCollapse<M,D,C,V,S,I>::Update_neighbors( vertex_pair_ptr aCollapsingPair ) 
{
  CGAL_TSMS_TRACE(3,"Updating cost of neighboring edges..." ) ;
  

  //
  // (A) Collect all pairs to update its cost: all those around each vertex adjacent to q  
  //
  
  typedef std::set<vertex_pair_ptr> vertex_pair_set ;
  
  vertex_pair_set lToUpdate ;
  
  
  // (A.1) Loop around all vertices adjacent to q
  in_edge_iterator eb1, ee1 ; 
  for ( tie(eb1,ee1) = in_edges(aCollapsingPair->q(),mSurface) ; eb1 != ee1 ; ++ eb1 )
  {
    edge_descriptor edge1 = *eb1 ;
    
    vertex_descriptor adj_q = source(edge1,mSurface);
    
    // (A.2) Loop around all edges incident on each adjacent vertex
    in_edge_iterator eb2, ee2 ; 
    for ( tie(eb2,ee2) = in_edges(adj_q,mSurface) ; eb2 != ee2 ; ++ eb2 )
    {
      edge_descriptor edge2 = *eb2 ;
      
      vertex_pair_ptr lPair = get_pair(edge2);
      
      if ( lPair->p() == aCollapsingPair->p() )
      {
        CGAL_TSMS_TRACE(4,"Replacing lPair->p() with Q" ) ;
        
        lPair->data().set(aCollapsingPair->q()
                         ,lPair          ->q()
                         ,aCollapsingPair->is_q_fixed() 
                         ,lPair          ->is_q_fixed()
                         ,source(edge2,mSurface) == aCollapsingPair->q() ? edge2 : opposite_edge(edge2,mSurface)
                         ,mSurface
                         );
      }
      else if ( lPair->q() == aCollapsingPair->p() )
      {
        CGAL_TSMS_TRACE(4,"Replacing lPair->q() with Q" ) ;
        
        lPair->data().set(lPair          ->p()
                         ,aCollapsingPair->q()
                         ,lPair          ->is_p_fixed()
                         ,aCollapsingPair->is_q_fixed() 
                         ,target(edge2,mSurface) == aCollapsingPair->q() ? edge2 : opposite_edge(edge2,mSurface)
                         ,mSurface
                         );
      }
      
      CGAL_TSMS_TRACE(4,"Inedge around V" << adj_q->ID << " E" << edge2->ID << " Opposite E" << edge2->opposite()->ID 
                     << " V" <<  edge2->opposite()->vertex()->ID << "->V" << edge2->vertex()->ID 
                     << "\nPair:" << *lPair
                     ) ;
                   
    
      // Only those pairs still in the PQ are update.
      // The mark is used because in the way we loop here the same pair is found many times.
      if ( lPair->is_in_PQ() )
      {
        if ( lToUpdate.find(lPair) == lToUpdate.end() )
        {
          CGAL_TSMS_TRACE(4,"Pair registered for updating.") ;
          lToUpdate.insert(lPair) ;
        }
      }  
    } 
  }  
  
  //
  // (B) Proceed to update the costs.
  //
  
  for ( typename vertex_pair_set::iterator it = lToUpdate.begin(), eit = lToUpdate.end() ; it != eit ; ++ it )
  {
    vertex_pair_ptr lPair = *it ;
    CGAL_TSMS_TRACE(3,"Updating cost of " << *lPair) ;
    
    // The cost of a pair can be recalculated by invalidating its cache and updating the PQ.
    // The PQ update will reposition the pair in the heap querying its cost(),
    // but since the cost was invalidated, it will be computed again 
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
    update_in_PQ(lPair);
  }
    
}


} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_IMPL_H //
// EOF //
 
