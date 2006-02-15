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
#ifndef CGAL_POLYGON_OFFSET_BUILDER_2_C
#define CGAL_POLYGON_OFFSET_BUILDER_2_C 1

CGAL_BEGIN_NAMESPACE

template<class Sls, class Gt, class Cont>
Polygon_offset_builder_2<Sls,Gt,Cont>::Polygon_offset_builder_2( Sls const& aSls, Traits const& aTraits )
  :
  mTraits(aTraits)
 ,mVisitedBisectors(aSls.size_of_halfedges())
{
  for ( Halfedge_const_handle lHE = aSls.halfedges_begin() ; lHE != aSls.halfedges_end() ; ++ lHE )
    if ( !lHE->is_bisector() )
      mBorders.push_back(lHE);

  ResetVisitedBisectorsMap();
}

template<class Sls, class Gt, class Cont>
typename Polygon_offset_builder_2<Sls,Gt,Cont>::Halfedge_const_handle
Polygon_offset_builder_2<Sls,Gt,Cont>::LocateHook( FT aTime, Halfedge_const_handle aBisector, bool aAbove )
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

template<class Sls, class Gt, class Cont>
typename Polygon_offset_builder_2<Sls,Gt,Cont>::Halfedge_const_handle
Polygon_offset_builder_2<Sls,Gt,Cont>::LocateSeed( FT aTime )
{
  Halfedge_const_handle rSeed ;

  for ( typename Halfedge_vector::const_iterator f = mBorders.begin()
       ; f != mBorders.end() && !handle_assigned(rSeed)
       ; ++ f
      )
    rSeed = LocateHook(aTime,(*f)->next(),true);

  return rSeed;
}

template<class Sls, class Gt, class Cont>
template<class OutputIterator>
OutputIterator Polygon_offset_builder_2<Sls,Gt,Cont>::TraceOffsetPolygon( FT aTime, Halfedge_const_handle aSeed, OutputIterator aOut )
{
  ContainerPtr lPoly( new Container() ) ;

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

template<class Sls, class Gt, class Cont>
void Polygon_offset_builder_2<Sls,Gt,Cont>::ResetVisitedBisectorsMap()
{
  std::fill(mVisitedBisectors.begin(),mVisitedBisectors.end(),0);
}

template<class Sls, class Gt, class Cont>
template<class OutputIterator>
OutputIterator Polygon_offset_builder_2<Sls,Gt,Cont>::construct_offset_polygons( FT aTime, OutputIterator aOut )
{
  ResetVisitedBisectorsMap();

  CGAL_SSBUILDER_TRACE("Constructing offset polygons for offset: " << aTime ) ;
  for ( Halfedge_const_handle lSeed = LocateSeed(aTime); handle_assigned(lSeed); lSeed = LocateSeed(aTime) )
    aOut = TraceOffsetPolygon(aTime,lSeed,aOut);

  return aOut ;
}

CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_OFFSET_BUILDER_2_C //
// EOF //
