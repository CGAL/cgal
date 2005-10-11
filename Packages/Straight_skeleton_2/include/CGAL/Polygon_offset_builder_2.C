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
// file          : include/CGAL/Polygon_offset_builder_2.c
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_POLYGON_OFFSET_BUILDER_2_C
#define CGAL_POLYGON_OFFSET_BUILDER_2_C 1

CGAL_BEGIN_NAMESPACE

template<class Sls, class Gt, class Poly>
Polygon_offset_builder_2<Sls,Gt,Poly>::Polygon_offset_builder_2( Sls const& aSls, Traits const& aTraits )
  :
  mTraits(aTraits)
 ,mVisitedBisectors(aSls.size_of_halfedges())
{
  for ( Halfedge_const_handle lHE = aSls.halfedges_begin() ; lHE != aSls.halfedges_end() ; ++ lHE )
    if ( !lHE->is_bisector() )
      mBorders.push_back(lHE);
}

template<class Sls, class Gt, class Poly>
typename Polygon_offset_builder_2<Sls,Gt,Poly>::Halfedge_const_handle
Polygon_offset_builder_2<Sls,Gt,Poly>::LocateHook( FT aTime, Halfedge_const_handle aBisector, bool aAbove )
{
  Halfedge_const_handle rHook ;

  while ( aBisector->is_bisector() )
  {
    Halfedge_const_handle lNext = aBisector->next();

    if ( !IsVisited(aBisector) )
    {
      if ( lNext->is_bisector() )
      {
        Comparison_result lC = Compare_offset_against_event_time(aTime,aBisector,lNext) ;
        if ( (aAbove && lC != SMALLER )  || (!aAbove && lC != LARGER) )
        {
          rHook = aBisector ;
          break ;
        }
      }
      else
      {
        if ( !aAbove )
          rHook = aBisector ;
      }

    }
    aBisector = lNext ;
  }

  return rHook;
}

template<class Sls, class Gt, class Poly>
typename Polygon_offset_builder_2<Sls,Gt,Poly>::Halfedge_const_handle
Polygon_offset_builder_2<Sls,Gt,Poly>::LocateSeed( FT aTime )
{
  Halfedge_const_handle rSeed ;

  for ( typename Halfedge_vector::const_iterator f = mBorders.begin()
       ; f != mBorders.end() && !handle_assigned(rSeed)
       ; ++ f
      )
    rSeed = LocateHook(aTime,(*f)->next(),true);

  return rSeed;
}

template<class Sls, class Gt, class Poly>
template<class OutputIterator>
OutputIterator Polygon_offset_builder_2<Sls,Gt,Poly>::TraceOffsetPolygon( FT aTime, Halfedge_const_handle aSeed, OutputIterator aOut )
{
  Polygon_2_Ptr lPoly( new Polygon_2() ) ;

  Halfedge_const_handle lHookL = aSeed ;

  while ( true )
  {
    //ACGEO_OFFSET_SHOW ( AutoACObject _( DrawSegment(lHookL.mBorderHE->segment(),Magenta,"OffsetAux") ) ) ;
    Visit(lHookL);

    Halfedge_const_handle lHookR = LocateHook(aTime,lHookL->next(),false) ;
    if ( handle_assigned(lHookR) )
    {
      Visit(lHookR);

      lPoly->push_back(Construct_offset_point(aTime,lHookR));

      Halfedge_const_handle lNextBisector = lHookR->opposite();

      if ( lNextBisector != aSeed && !IsVisited(lNextBisector) )
      {
        lHookL = lNextBisector;
        continue;
      }
    }
    break ;
  }

  CGAL_SSBUILDER_TRACE("Offset polygon of " << lPoly->size() << " vertices traced." ) ;

  if ( lPoly->size() >= 3 )
    *aOut++ = lPoly ;

  return aOut ;
}

template<class Sls, class Gt, class Poly>
template<class OutputIterator>
OutputIterator Polygon_offset_builder_2<Sls,Gt,Poly>::Create( FT aTime, OutputIterator aOut )
{
  mVisitedBisectors.clear();
  CGAL_SSBUILDER_TRACE("Constructing offset polygons for offset: " << aTime ) ;
  for ( Halfedge_const_handle lSeed = LocateSeed(aTime); handle_assigned(lSeed); lSeed = LocateSeed(aTime) )
    aOut = TraceOffsetPolygon(aTime,lSeed,aOut);

  return aOut ;
}

CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_OFFSET_BUILDER_2_C //
// EOF //
