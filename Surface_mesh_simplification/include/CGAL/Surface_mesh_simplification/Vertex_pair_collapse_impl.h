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

template<class M,class D,class C,class V,class S>
VertexPairCollapse<M,D,C,V,S>::VertexPairCollapse( TSM&                           aSurface
                                                 , GetCollapseData const&         aGet_collapse_data
                                                 , ParamsToGetCollapseData const* aParamsToGetCollapseData
                                                 , GetCost const&                 aGet_cost
                                                 , GetNewVertexPoint const&       aGet_new_vertex_point
                                                 , ShouldStop const&              aShould_stop
                                                 , bool                           aIncludeNonEdgePairs 
                                                 )
  : 
   mSurface (aSurface)
  ,mParamsToGetCollapseData(aParamsToGetCollapseData)
   
  ,Get_collapse_data   (aGet_collapse_data)
  ,Get_cost            (aGet_cost)
  ,Get_new_vertex_point(aGet_new_vertex_point)
  ,Should_stop         (aShould_stop) 
  
  ,mIncludeNonEdgePairs(aIncludeNonEdgePairs)
{
  CGAL_TSMS_TRACE(0,"VertexPairCollapse of TSM with " << num_edges(aSurface)/2 << " edges" );
}

template<class M,class D,class C,class V,class S>
int VertexPairCollapse<M,D,C,V,S>::run()
{
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
    mInitialPairCount = mCurrentPairCount = num_edges(mSurface) / 2 ;
    
    CGAL_TSMS_TRACE(0,"A thetrahedron cannot be simplified any further.");
  }
  
  return (int)(mInitialPairCount - mCurrentPairCount) ;
}

template<class M,class D,class C,class V,class S>
void VertexPairCollapse<M,D,C,V,S>::Collect()
{
  CGAL_TSMS_TRACE(0,"Collecting vertex-pairs...");

  // Loop over all the edges in the surface putting the accepted vertex-pairs in the PQ
  
  Equal_3 equal_points = Kernel().equal_3_object();
    
  size_type lID = 0 ;
  
  undirected_edge_iterator eb, ee ;
  for ( tie(eb,ee) = undirected_edges(mSurface); eb!=ee; ++eb )
  {
    edge_descriptor edge = *eb ;
    
    vertex_descriptor s = source(edge,mSurface);
    vertex_descriptor t = target(edge,mSurface);
    
    Point_3 sp = get_point(s) ;
    Point_3 tp = get_point(t) ;
    
    if ( ! equal_points(sp,tp) )
    {
      bool is_s_fixed, is_t_fixed ;
      GetVerticesIsFixedFlags(s,t,edge,is_s_fixed,is_t_fixed);
 
      // The collapse always keeps the vertex 't' and removes vertex 's'.
      // If the 's' is fixed (but not t) then swap the pair to keep (the original) 's' correctly fixed.
      if ( is_s_fixed && !is_t_fixed )
      {
        using std::swap ;
        swap(s,t);
        swap(is_s_fixed,is_t_fixed);
        edge = opposite_edge(edge,mSurface);
      }
      
      Collapse_data_ptr lData = Get_collapse_data(s,t,is_s_fixed,is_t_fixed,edge,mSurface,mParamsToGetCollapseData) ;
      vertex_pair_ptr lPair( new vertex_pair(lID++,lData,Get_cost) ) ; 
      set_pair(edge,lPair);
      insert_in_PQ(lPair);
      CGAL_TSMS_TRACE(3, (lPair->cost(),*lPair) << " accepted." );
      CGAL_TSMS_AUDIT(s,t,edge,lPair->cost(),Get_new_vertex_point(*lPair->data()));
    }
  }
  

  // TODO: Collect non-edge pairs if requested 
  
  mInitialPairCount = mCurrentPairCount = lID;
  
  CGAL_TSMS_TRACE(0,"Initial pair count: " << mInitialPairCount ) ;
}

template<class M,class D,class C,class V,class S>
void VertexPairCollapse<M,D,C,V,S>::Loop()
{
  CGAL_TSMS_TRACE(0,"Removing pairs...") ;

  //
  // Pops and processes each vertex-pair from the PQ
  //
  while ( vertex_pair_ptr lPair = pop_from_PQ() )
  {
    CGAL_TSMS_TRACE(3,"Poped " << *lPair ) ;
    
    if ( lPair->cost() != none ) 
    {
      if ( num_vertices(mSurface) <= 4 )
      {
        CGAL_TSMS_TRACE(0,"Thetrahedron reached.");
        break ;
      }
      
      if ( Should_stop(*lPair->cost(),*lPair->data(),mInitialPairCount,mCurrentPairCount) )
      {
        CGAL_TSMS_TRACE(0,"Stop condition reached with InitialCount=" << mInitialPairCount
                       << " CurrentPairCount=" << mCurrentPairCount
                       << " CurrentPair: " << *lPair
                       );
        break ;
      }

      Collapse(lPair);
    }
  }
}

// Some edges are NOT collapsable: doing so would break the topological consistency of the mesh.
// This function returns true if a edge 'p->q' can be collapsed.
//
template<class M,class D,class C,class V,class S>
bool VertexPairCollapse<M,D,C,V,S>::Is_collapsable( vertex_descriptor const& p, vertex_descriptor const& q, edge_descriptor const& p_q )
{
  bool rR = true ;
  
  edge_descriptor q_p = opposite_edge(p_q,mSurface);
  
  edge_descriptor p_t = next_edge_ccw(p_q,mSurface);      
  edge_descriptor p_b = next_edge_cw (p_q,mSurface);      
  edge_descriptor q_t = next_edge_cw (q_p,mSurface);      
  edge_descriptor q_b = next_edge_ccw(q_p,mSurface);      

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
    else rR = false ;
  }   
  else rR = false ;
  
  return rR ;
}

template<class M,class D,class C,class V,class S>
void VertexPairCollapse<M,D,C,V,S>::Collapse( vertex_pair_ptr const& aPair )
{
  CGAL_TSMS_TRACE(1,"Collapsig " << *aPair ) ;
  
  vertex_descriptor lP = aPair->p();
  vertex_descriptor lQ = aPair->q();
  
  CGAL_assertion( lP != lQ );
  
  // This external function is allowed to return an absent point if there is no way to place the vertex
  // satisfying its constrians. In that case the vertex-pair is simply not removed.
  Optional_vertex_point_type lNewVertexPoint = Get_new_vertex_point(*aPair->data());
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
                   << " V" <<  lEdgePQ->opposite()->vertex()->ID << "->V" << lEdgePQ->vertex()->ID ) ;
    CGAL_TSMS_TRACE(3,"EdgePT E" << lEdgePT->ID << " Opposite EdgePT E" << lEdgePT->opposite()->ID 
                   << " V" <<  lEdgePT->opposite()->vertex()->ID << "->V" << lEdgePT->vertex()->ID ) ;

    vertex_pair_vector lUpdated ;
    
    // Since vertex P will be removed during the collapse, all cached pairs linked to 'P' (either from or to) 
    // are updated to link to 'Q' instead.
    out_edge_iterator eb, ee ; 
    for ( boost::tie(eb,ee) = boost::out_edges(lP,mSurface) ; eb != ee ; ++ eb )
    {
      edge_descriptor outedge = *eb ;
      CGAL_TSMS_TRACE(4,"Outedge around V" << lP->ID << " E" << outedge->ID << " Opposite E" << outedge->opposite()->ID 
                   << " V" <<  outedge->opposite()->vertex()->ID << "->V" << outedge->vertex()->ID ) ;
      
      if ( outedge != lEdgePQ && outedge != lEdgePT )
      {
        vertex_pair_ptr lPair = get_pair(outedge);
        CGAL_TSMS_TRACE(4,"Updating vertex P in " << *lPair) ;
        
        CGAL_assertion( lPair->p() == lP || lPair->q() == lP );
      
        bool is_v0_fixed, is_v1_fixed ;
        vertex_descriptor v0, v1 ;
        edge_descriptor   edge ;
        if ( lPair->p() == lP )
        {
          v0 = lQ ;
          v1 = lPair->q() ;
          is_v0_fixed = aPair->is_q_fixed() ;
          is_v1_fixed = lPair->is_q_fixed() ;
          edge = outedge ;
        }
        else
        {
          v0 = lPair->p() ;
          v1 = lQ ;
          is_v0_fixed = lPair->is_p_fixed() ;
          is_v1_fixed = aPair->is_q_fixed() ;
          edge = opposite_edge(outedge,mSurface);
        }
      
        if ( !Is_collapsable(v0,v1,edge) )
          is_v0_fixed = is_v1_fixed = true ;
        
        Collapse_data_ptr lNewData = Get_collapse_data(v0,v1,is_v0_fixed,is_v1_fixed,edge,mSurface,mParamsToGetCollapseData);
        lPair->reset_data(lNewData);
        lPair->mark() = 1 ;
        lUpdated.push_back(lPair);
        
        CGAL_TSMS_TRACE(4,"...after update: " << *lPair ) ;
      }
    }
    
    // The collapse will remove the following edges from the surface so the corresponding pairs won't be valid anymore.
    
    vertex_pair_ptr lPairPT = get_pair(lEdgePT) ;
    vertex_pair_ptr lPairQB = get_pair(lEdgeQB) ;
    
    if ( lPairPT->is_in_PQ() )
    {
      CGAL_TSMS_TRACE(2,"Removing from PQ VP" << lPairPT->id() ) ;
      remove_from_PQ(lPairPT) ;
    }
    
    if ( lPairQB->is_in_PQ() )
    {
      CGAL_TSMS_TRACE(2,"Removing from PQ VP" << lPairQB->id() ) ;
      remove_from_PQ(lPairQB) ;
    }
    
    CGAL_TSMS_TRACE(2,"Removing from surface V" << lP->ID 
                   << " E" << lEdgePQ->ID
                   << " E" << lEdgePT->ID
                   << " E" << lEdgeQB->ID
                   );

                         
    Collapse_triangulation_edge(lEdgePQ,lEdgePT,lEdgeQB,mSurface);
    
    // Reset the point of placement of Q (the vertex that "replaces" the collapsed edge)
    vertex_point_t vertex_point ;
    put(vertex_point,mSurface,lQ,*lNewVertexPoint) ;

    // Updates the cost of all pairs in the PQ
    Update_neighbors(lQ);
            
    for ( vertex_pair_vector_iterator it = lUpdated.begin(), eit = lUpdated.end() ; it != eit ; ++ it )
      (*it)->mark() = 0 ;
      
    mCurrentPairCount -= 3 ;
  }
  else
  {
    CGAL_TSMS_TRACE(0,"Unable to calculate new vertex point. Pair " << aPair << " discarded without being removed" ) ;
  }
}


template<class M,class D,class C,class V,class S>
void VertexPairCollapse<M,D,C,V,S>::Update_neighbors( vertex_descriptor const& v ) 
{
  CGAL_TSMS_TRACE(3,"Updating cost of neighboring edges..." ) ;
  

  //
  // (A) Collect all pairs to update its cost: all those around each vertex adjacent to v  
  //
  
  vertex_pair_vector lToUpdate ;
  
  // (A.1) Loop around all vertices adjacent to v
  in_edge_iterator eb1, ee1 ; 
  for ( tie(eb1,ee1) = in_edges(v,mSurface) ; eb1 != ee1 ; ++ eb1 )
  {
    edge_descriptor edge1 = *eb1 ;
    
    vertex_pair_ptr lPair1 = get_pair(edge1) ;
    
    CGAL_TSMS_TRACE(4,"Inedge around V" << v->ID << " E" << edge1->ID << " Opposite E" << edge1->opposite()->ID 
                   << " V" <<  edge1->opposite()->vertex()->ID << "->V" << edge1->vertex()->ID
                   << "\n" << *lPair1
                   ) ;
                   
    
    vertex_descriptor adj_v = source(edge1,mSurface);
    
    // (A.2) Loop around all edges incident on each adjacent vertex
    in_edge_iterator eb2, ee2 ; 
    for ( tie(eb2,ee2) = in_edges(adj_v,mSurface) ; eb2 != ee2 ; ++ eb2 )
    {
      edge_descriptor edge2 = *eb2 ;
      
      vertex_pair_ptr lPair2 = get_pair(edge2);
      
      CGAL_TSMS_TRACE(4,"Inedge around V" << adj_v->ID << " E" << edge2->ID << " Opposite E" << edge2->opposite()->ID 
                     << " V" <<  edge2->opposite()->vertex()->ID << "->V" << edge2->vertex()->ID 
                     << "\n" << *lPair2
                     ) ;
                   
    
      // Only those pairs still in the PQ are update.
      // The mark is used because in the way we loop here the same pair is found many times.
      if ( lPair2->is_in_PQ() && lPair2->mark() == 0 )
      {
        CGAL_TSMS_TRACE(4,"Pair registered for updating.") ;
        lPair2->mark() = 1 ;
        lToUpdate.push_back(lPair2);
      }  
    } 
  }  
  
  //
  // (B) Proceed to update the costs.
  //
  
  for ( vertex_pair_vector_iterator it = lToUpdate.begin(), eit = lToUpdate.end() ; it != eit ; ++ it )
  {
    vertex_pair_ptr lPair = *it ;
    CGAL_TSMS_TRACE(3,"Updating cost of " << *lPair) ;
    
    // The cost of a pair can be recalculated by invalidating its cache and updating the PQ.
    // The PQ update will reposition the pair in the heap querying its cost(),
    // but since the cost was invalidated, it will be computed again 
    Collapse_data_ptr lNewData = Get_collapse_data(lPair->p()
                                                  ,lPair->q()
                                                  ,lPair->is_p_fixed()
                                                  ,lPair->is_q_fixed()
                                                  ,lPair->edge()
                                                  ,lPair->surface()
                                                  ,mParamsToGetCollapseData
                                                  ) ;
    lPair->reset_data(lNewData);
    update_in_PQ(lPair);
    lPair->mark() = 0 ;
  }
    
}

} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_IMPL_H //
// EOF //
 
