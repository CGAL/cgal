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
// $URL: $
// $Id: $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_C
#define CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_C

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{

template<class TSM,class SM,class CM,class VP,class SP>
VertexPairCollapse<TSM,SM,CM,VP,SP>::VertexPairCollapse( TSM&                   aSurface
                                                       , SelectMap       const& aSelectMap
                                                       , CostMap         const& aCostMap
                                                       , VertexPlacement const& aVertexPlacement
                                                       , StopPred        const& aStopPred 
                                                       , bool                   aIncludeNonEdgePairs 
                                                       )
  : 
   mSurface (aSurface)
   
  ,Select_map(aSelectMap)
  ,Cost_map(aCostMap)
  
  ,construct_new_vertex_point(aVertexPlacement)
  ,stop_simplification(aStopPred) 
  
  ,mIncludeNonEdgePairs(aIncludeNonEdgePairs)
{
  CGAL_TSMS_TRACE(0,"VertexPairCollapse of TSM with " << boost::num_edges(aSurface) << " edges" );
}

template<class TSM,class SM,class CM,class VP,class SP>
int VertexPairCollapse<TSM,SM,CM,VP,SP>::run()
{
  Collect(); 
  Loop();
  
  CGAL_TSMS_TRACE(0,"Finished: " << (mInitialPairCount - mCurrentPairCount) << " pairs removed." ) ;
  
  return (int)(mInitialPairCount - mCurrentPairCount) ;
}

template<class TSM,class SM,class CM,class VP,class SP>
void VertexPairCollapse<TSM,SM,CM,VP,SP>::Collect()
{
  CGAL_TSMS_TRACE(0,"Collecting vertex-pairs...");

  // Loop over all the edges in the surface putting the accepted vertex-pairs in the PQ
  
  boost::vertex_point_t vertex_point ;
  typedef typename Kernel_traits<Point_3>::Kernel Kernel ;
  Kernel kernel ;
  typename Kernel::Equal_3 equal_points = kernel.equal_3_object();
    
  size_type lID = 0 ;
  
  edge_iterator eb, ee ;
  for ( boost::tie(eb,ee) = boost::edges(mSurface); eb!=ee; ++eb )
  {
    edge_descriptor edge = *eb ;
    
    vertex_descriptor s = boost::source(edge,mSurface);
    vertex_descriptor t = boost::target(edge,mSurface);
    
    Point_3 sp = boost::get(vertex_point,mSurface,s) ;
    Point_3 tp = boost::get(vertex_point,mSurface,t) ;
    
    if ( ! equal_points(sp,tp) )
    {
      if ( boost::get(Select_map, boost::make_tuple(s,t,true,boost::addressof(mSurface))) )
      {
        vertex_pair_ptr lPair( new vertex_pair(lID++,s,t,edge,boost::addressof(mSurface),boost::addressof(Cost_map)) ) ; 
        set_pair(edge,lPair);
        insert_in_PQ(lPair);
        CGAL_TSMS_TRACE(3, *lPair << " accepted." );
      }    
    }
  }
  

  // TODO: Collect non-edge pairs if requested 
  
  mInitialPairCount = mCurrentPairCount = lID;
  
  CGAL_TSMS_TRACE(0,"Initial pair count: " << mInitialPairCount ) ;
}

template<class TSM,class SM,class CM,class VP,class SP>
void VertexPairCollapse<TSM,SM,CM,VP,SP>::Loop()
{
  CGAL_TSMS_TRACE(0,"Removing pairs...") ;

  //
  // Pops and processes each vertex-pair from the PQ
  //
  while ( vertex_pair_ptr lPair = pop_from_PQ() )
  {
    CGAL_TSMS_TRACE(3,"Poped " << *lPair ) ;
    
    // Pairs in the queue might be "fixed", that is, marked as uncollapsable, or their cost might be undefined.
    if ( !lPair->is_fixed() && lPair->cost() != boost::none ) 
    {
      if ( boost::num_vertices(mSurface) <= 4 )
      {
        CGAL_TSMS_TRACE(0,"Thetrahedron reached.");
        break ;
      }
      
      if ( stop_simplification(*lPair->cost(),lPair->p(),lPair->q(),lPair->is_edge(),mInitialPairCount,mCurrentPairCount,mSurface) )
      {
        CGAL_TSMS_TRACE(0,"Stop condition reached with InitialCount=" << mInitialPairCount
                       << " CurrentPairCount=" << mCurrentPairCount
                       << " CurrentPair: " << *lPair
                       );
        break ;
      }

     // Proceeds only if the pair is topolofically collapsable (collapsing it results in a non-degenerate triangulation)
     if ( Is_collapsable(lPair) )        
        Collapse(lPair);
    }
  }
}

template<class TSM,class SM,class CM,class VP,class SP>
bool VertexPairCollapse<TSM,SM,CM,VP,SP>::Is_collapsable( vertex_pair_ptr const& aPair )
{
  bool rR = true ;
  
  edge_descriptor p_q = aPair->edge();
  edge_descriptor q_p = boost::opposite_edge(p_q,mSurface);
  
  edge_descriptor p_t = boost::next_edge_ccw(p_q,mSurface);      
  edge_descriptor p_b = boost::next_edge_cw (p_q,mSurface);      
  edge_descriptor q_t = boost::next_edge_cw (q_p,mSurface);      
  edge_descriptor q_b = boost::next_edge_ccw(q_p,mSurface);      

  // degree(p) and degree(q) >= 3
  if (    boost::target(p_t,mSurface) != aPair->q()
       && boost::target(p_b,mSurface) != aPair->q()
       && boost::target(q_t,mSurface) != aPair->p()
       && boost::target(q_b,mSurface) != aPair->p()
     )
  {
    // link('p') .intersection. link('q') == link('p_q') (that is, exactly {'t','b'})
    if (    boost::target(p_t,mSurface) == boost::target(q_t,mSurface) 
         && boost::target(p_b,mSurface) == boost::target(q_b,mSurface) 
       )
    {
      edge_descriptor p_tn = p_t ;
      edge_descriptor p_bn = p_b ;
      edge_descriptor q_tn = q_t ;
      edge_descriptor q_bn = q_b ;
      
      do
      {
        p_tn = boost::next_edge_ccw(p_tn,mSurface);      
        p_bn = boost::next_edge_cw (p_bn,mSurface);      
        q_tn = boost::next_edge_cw (q_tn,mSurface);      
        q_bn = boost::next_edge_ccw(q_bn,mSurface);      
        
        if (    boost::target(p_tn,mSurface) == boost::target(q_tn,mSurface)
             || boost::target(p_bn,mSurface) == boost::target(q_bn,mSurface)
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

template<class TSM,class SM,class CM,class VP,class SP>
void VertexPairCollapse<TSM,SM,CM,VP,SP>::Collapse( vertex_pair_ptr const& aPair )
{
  CGAL_TSMS_TRACE(1,"Collapsig " << *aPair ) ;
  
  vertex_descriptor lP = aPair->p();
  vertex_descriptor lQ = aPair->q();
  
  // This external function is allowed to return an absent point if there is no way to place the vertex
  // satisfying its constrians. In that case the vertex-pair is simply not removed.
  optional_vertex_point_type lNewVertexPoint = construct_new_vertex_point(lP,lQ,aPair->is_edge(),mSurface);
  if ( lNewVertexPoint )
  {
    CGAL_TSMS_TRACE(2,"New vertex point: (" << lNewVertexPoint->x() << "," << lNewVertexPoint->y() << "," << lNewVertexPoint->z() << ")");
        
    // The actual collapse of edge PQ merges the top and bottom facets with its left and right adjacents resp, then
    // joins P and Q.
    
    edge_descriptor lEdgePQ = aPair->edge();
    
    edge_descriptor lEdgeQP = boost::opposite_edge(lEdgePQ,mSurface);
    
    edge_descriptor lEdgePT = boost::next_edge_ccw(lEdgePQ,mSurface);
    
    edge_descriptor lEdgeQB = boost::next_edge_ccw(lEdgeQP,mSurface);
    
    CGAL_TSMS_TRACE(3,"EdgePQ E" << lEdgePQ->ID << " Opposite EdgePQ E" << lEdgePQ->opposite()->ID 
                   << " V" <<  lEdgePQ->opposite()->vertex()->ID << "->V" << lEdgePQ->vertex()->ID ) ;
    CGAL_TSMS_TRACE(3,"EdgePT E" << lEdgePT->ID << " Opposite EdgePT E" << lEdgePT->opposite()->ID 
                   << " V" <<  lEdgePT->opposite()->vertex()->ID << "->V" << lEdgePT->vertex()->ID ) ;

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
        lPair->update_vertex(lP,lQ) ;
        CGAL_TSMS_TRACE(4,"...after update: " << *lPair ) ;
      }
    }
    
    vertex_pair_ptr lPairPT = get_pair(lEdgePT) ;
    vertex_pair_ptr lPairQB = get_pair(lEdgeQB) ;
  
    // The collapse will remove these edges from the surface so the corresponding pairs won't be valid anymore.
    
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
    
    // This operator IS NOT passed as a policy. It is a traits. Users only need to specialize it for
    // the particular surface type.
    collapse_triangulation_edge(lEdgePQ,lEdgePT,lEdgeQB,mSurface);
    
    // Reset the point of placement of Q (the vertex that "replaces" the collapsed edge)
    boost::vertex_point_t vertex_point ;
    boost::put(vertex_point,mSurface,lQ,*lNewVertexPoint) ;

    // Updates the cost of all pairs in the PQ
    Update_neighbors(lQ);
            
    mCurrentPairCount -= 3 ;
  }
  else
  {
    CGAL_TSMS_TRACE(0,"Unable to calculate new vertex point. Pair " << aPair << " discarded without being removed" ) ;
  }
}


template<class TSM,class SM,class CM,class VP,class SP>
void VertexPairCollapse<TSM,SM,CM,VP,SP>::Update_neighbors( vertex_descriptor const& v ) 
{
  CGAL_TSMS_TRACE(3,"Updating cost of neighboring edges..." ) ;
  

  //
  // (A) Collect all pairs to update its cost: all those around each vertex adjacent to v  
  //
  
  vertex_pair_vector lToUpdate ;
  
  // (A.1) Loop around all vertices adjacent to v
  in_edge_iterator eb1, ee1 ; 
  for ( boost::tie(eb1,ee1) = boost::in_edges(v,mSurface) ; eb1 != ee1 ; ++ eb1 )
  {
    edge_descriptor edge1 = *eb1 ;
    
    CGAL_TSMS_TRACE(4,"Inedge around V" << v->ID << " E" << edge1->ID << " Opposite E" << edge1->opposite()->ID 
                   << " V" <<  edge1->opposite()->vertex()->ID << "->V" << edge1->vertex()->ID ) ;
                   
    vertex_pair_ptr lPair1 = get_pair(edge1) ;
    
    // This is required to satisfy the transitive link_condition.
    // That is, the edges around the replacement vertex 'v' cannot be collapsed again.
    //lPair1->is_fixed() = true ;
   
    vertex_descriptor adj_v = boost::source(edge1,mSurface);
    
    // (A.2) Loop around all edges incident on each adjacent vertex
    in_edge_iterator eb2, ee2 ; 
    for ( boost::tie(eb2,ee2) = boost::in_edges(adj_v,mSurface) ; eb2 != ee2 ; ++ eb2 )
    {
      edge_descriptor edge2 = *eb2 ;
      
      CGAL_TSMS_TRACE(4,"Inedge around V" << adj_v->ID << " E" << edge2->ID << " Opposite E" << edge2->opposite()->ID 
                   << " V" <<  edge2->opposite()->vertex()->ID << "->V" << edge2->vertex()->ID ) ;
                   
      vertex_pair_ptr lPair2 = get_pair(edge2);
    
      // Only those pairs still in the PQ are update.
      // The mark is used because in the way we loop here the same pair is found many times.
      if ( lPair2->is_in_PQ() && lPair2->mark() == 0 )
      {
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
    lPair->invalidate_cost();
    update_in_PQ(lPair);
    lPair->mark() = 0 ;
  }
    
}

} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_C //
// EOF //
 
