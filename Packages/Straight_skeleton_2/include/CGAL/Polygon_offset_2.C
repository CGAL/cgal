// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Polygon_offset_2.c
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_POLYGON_OFFSET_2_C
#define CGAL_POLYGON_OFFSET_2_C 1


CGAL_BEGIN_NAMESPACE

template<class Graph, class Traits>
Polygon_offset_2<Graph,Traits>::Polygon_offset_2( Graph const& aGraph, Traits const& aTraits )
  :
  mTraits(aTraits)
 ,mVisitedBisectors(aGraph.size_of_halfedges())   
{
  for ( Halfedge_const_handle lHE = aGraph.halfedges_begin() ; lHE != aGraph.halfedges_end() ; ++ lHE )
    if ( !lHE->is_bisector() )
      mBorders.push_back(lHE);
}

template<class Graph, class Traits>
typename Polygon_offset_2<Graph,Traits>::Halfedge_const_handle
Polygon_offset_2<Graph,Traits>::LocateHook( FT aTime, Halfedge_const_handle aBisector )
{
  Halfedge_const_handle rHook ;
  
  while ( aBisector->is_bisector() )
  {
    if ( !IsVisited(aBisector) )
    {
      //ACGEO_OFFSET_SHOW( AutoACObject _( DrawSegment(aBisectHE->segment(),Cyan,"OffsetAux") ) ) ;
      if ( IsTimeWithinBisector(aTime,aBisector) )
      {
        rHook = aBisector ;
        break ;
      }
    }    
    aBisector = aBisector->next();    
  }
         
  return rHook;
}

template<class Graph, class Traits>
typename Polygon_offset_2<Graph,Traits>::Halfedge_const_handle
Polygon_offset_2<Graph,Traits>::LocateSeed( FT aTime )
{
  Halfedge_const_handle rSeed ;
  
  for ( typename Halfedge_vector::const_iterator f = mBorders.begin()
       ; f != mBorders.end() && !handle_assigned(rSeed)
       ; ++ f
      )
    rSeed = LocateHook(aTime,(*f)->next());  
       
  return rSeed;
}

template<class Graph, class Traits>
template<class OutputIterator>
OutputIterator Polygon_offset_2<Graph,Traits>::TraceOffsetPolygon( FT                    aTime
                                                                 , Halfedge_const_handle aSeed
								 , OutputIterator        aOut 
								 )
{
  
  Polygon_2_Ptr lPoly = mTraits.Construct_polygon_2_ptr();
 
  Halfedge_const_handle lHookL = aSeed ;
  
  while ( true )
  {
    //ACGEO_OFFSET_SHOW ( AutoACObject _( DrawSegment(lHookL.mBorderHE->segment(),Magenta,"OffsetAux") ) ) ;
    Visit(lHookL);
    
    Halfedge_const_handle lHookR = LocateHook(aTime,lHookL->next()) ;
    if ( handle_assigned(lHookR) )
    {
      Visit(lHookR); 
      
      CGAL_assertion(handle_assigned(lHookR->opposite()));
      //ACGEO_OFFSET_SHOW ( DrawSegment(lHookL.mP,lHookR->mP,Blue,"OffsetPolygon") ) ;
      
      Halfedge_const_handle lBorder = lHookR->defining_border();
      
      CGAL_assertion(handle_assigned(lBorder));
      CGAL_assertion(handle_assigned(lBorder->opposite()));
      
      Optional_point_2 lP = mTraits.Construct_offset_point_2(aTime
                                                            ,lBorder->opposite()->vertex()->point() 
                                                            ,lBorder            ->vertex()->point() 
                                                            ,lHookR->opposite() ->vertex()->point()
                                                            ,lHookR             ->vertex()->point()
                                                            ) ;
      if ( lP )
        lPoly->push_back(*lP);
      
      Halfedge_const_handle lNextBisector = lHookR->opposite();
      CGAL_assertion(handle_assigned(lNextBisector));
      
      if ( lNextBisector != aSeed && !IsVisited(lNextBisector) )
      {
        lHookL = lNextBisector; 
        continue;
      }
    }
    break ;
  } 
  
  *aOut++ = lPoly ;
    
  return aOut ;  
}

template<class Graph, class Traits>
template<class OutputIterator>
OutputIterator Polygon_offset_2<Graph,Traits>::Create( FT aTime, OutputIterator aOut )
{
  for ( Halfedge_const_handle lSeed = LocateSeed(aTime)  
      ; handle_assigned(lSeed)
      ; lSeed = LocateSeed(aTime) 
      )
    aOut = TraceOffsetPolygon(aTime,lSeed,aOut);
    
  return aOut ;  
}

CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_OFFSET_2_C //
// EOF //
