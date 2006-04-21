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
  
  ,mPQ( boost::num_edges(aSurface) )
{
  CGAL_TSMS_TRACE(0,"VertexPairCollapse of TSM with " << boost::num_edges(aSurface) << " edges" );
}

template<class TSM,class SM,class CM,class VP,class SP>
std::size_t VertexPairCollapse<TSM,SM,CM,VP,SP>::run()
{
  Collect(); 
  Loop();
  
  return mInitialPairCount - mCurrentPairCount ;
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
        edge_descriptor opposite = boost::opposite_edge(edge,mSurface);
        
        vertex_pair_ptr lPair( new vertex_pair(lID++,s,t,edge,boost::addressof(mSurface),boost::addressof(Cost_map)) ) ; 
        CGAL_TSMS_TRACE(2, *lPair << " accepted." );
        mEdgesToPairsMap[edge    ]=lPair;
        mEdgesToPairsMap[opposite]=lPair;
        mPQ.push(lPair);
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
  while ( !mPQ.empty() )
  {
    vertex_pair_ptr lPair = mPQ.top();
    mPQ.pop();
    
    // The cost of a pair in the queue might be left undefined after a previous collapse
    if ( lPair->cost() != boost::none ) 
    { 
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
    else
    {
      CGAL_TSMS_TRACE(0,lPair << " has absent cost now. Discarded without being removed.");
    }
  }
}

template<class TSM,class SM,class CM,class VP,class SP>
bool VertexPairCollapse<TSM,SM,CM,VP,SP>::Is_collapsable( vertex_pair_ptr const& aPair )
{
  return    boost::degree(aPair->p(),mSurface) > 3
         && boost::degree(aPair->q(),mSurface) > 3 ;
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
    CGAL_TSMS_TRACE(1,"New vertex point: (" << lNewVertexPoint->x() << "," << lNewVertexPoint->y() << "," << lNewVertexPoint->z() << ")");
        
    // The actual collapse of edge PQ merges the top and bottom facets with its left and right adjacents resp, then
    // joins P and Q.
    
    edge_descriptor lEdgePQ = aPair->edge();
    
    edge_descriptor lEdgePT = boost::next_out_edge_ccw(lEdgePQ,mSurface);
    
    edge_descriptor lEdgeQB = boost::next_out_edge_ccw(boost::opposite_edge(lEdgePQ,mSurface),mSurface);
    
    CGAL_TSMS_TRACE(1,"Removing:\n  " << *get_pair(lEdgePQ) << "\n  " << *get_pair(lEdgePT) << "\n  " << *get_pair(lEdgeQB) ) ;
    
    // The collapse will remove the edges from the surface so the corresponding pairs won't be valid anymore.
    mPQ.remove( get_pair(lEdgePQ) ) ;
    mPQ.remove( get_pair(lEdgePT) ) ;
    mPQ.remove( get_pair(lEdgeQB) ) ;

    // Since vertex P will be removed during the collapse, all cached pairs linked to 'P' (either from or to) 
    // are updated to link to 'Q' instead.
    in_edge_iterator eb, ee ; 
    for ( boost::tie(eb,ee) = boost::in_edges(lP,mSurface) ; eb != ee ; ++ eb )
      get_pair(*eb)->update_vertex(lP,lQ) ;
    
    // This operator IS NOT passed as a policy. It is a traits. Users only need to specialize it for
    // the particular surface type.
    collapse_triangulation_edge(lEdgePQ,lEdgePT,lEdgeQB,mSurface);
    
    // Reset the point of placement of Q (the vertex that "replaces" the collapsed edge)
    boost::vertex_point_t vertex_point ;
    boost::put(vertex_point,mSurface,lQ,*lNewVertexPoint) ;

    // Updates the cost of all pairs in the PQ
    Update_neighbors(lQ);
            
    -- mCurrentPairCount ;
  }
  else
  {
    CGAL_TSMS_TRACE(0,"Unable to calculate new vertex point. Pair " << aPair << " discarded without being removed" ) ;
  }
}


template<class TSM,class SM,class CM,class VP,class SP>
void VertexPairCollapse<TSM,SM,CM,VP,SP>::Update_neighbors( vertex_descriptor const& v ) 
{
  CGAL_TSMS_TRACE(1,"Updating cost of neighboring edges..." ) ;
  
  //
  // The cost of each edge around each vertex adjacent to new vertex is updated.
  //
  
  // (1) Loop around all vertices adjacent to v
  in_edge_iterator eb1, ee1 ; 
  for ( boost::tie(eb1,ee1) = boost::in_edges(v,mSurface) ; eb1 != ee1 ; ++ eb1 )
  {
    vertex_descriptor adj_v = boost::source(*eb1,mSurface);
    
    // (2) Loop around all edges incident on each adjacent vertex
    in_edge_iterator eb2, ee2 ; 
    for ( boost::tie(eb2,ee2) = boost::in_edges(adj_v,mSurface) ; eb2 != ee2 ; ++ eb2 )
    {
      vertex_pair_ptr lPair = get_pair(*eb2);
      
      CGAL_TSMS_TRACE(4,"Updating cost of " << *lPair) ;
      
      // The cost of a pair can be recalculated by invalidating its cache and updating the PQ.
      // The PQ update will reposition the pair in the heap querying its cost(),
      // but since the cost was invalidated, it will be computed again 
      // (that's why the pairs hold a pointer to the CostMap)
      lPair->invalidate_cost();
      mPQ.update(lPair);
    }  
  }  
}

} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_C //
// EOF //
 
