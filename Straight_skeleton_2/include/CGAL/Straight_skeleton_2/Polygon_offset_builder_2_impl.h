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
#ifndef CGAL_POLYGON_OFFSET_BUILDER_2_IMPL_H
#define CGAL_POLYGON_OFFSET_BUILDER_2_IMPL_H 1

CGAL_BEGIN_NAMESPACE

template<class Ss, class Gt, class Cont>
Polygon_offset_builder_2<Ss,Gt,Cont>::Polygon_offset_builder_2( Ss const& aSs, Traits const& aTraits )
  :
  mTraits(aTraits)
{

  int lMaxID = -1 ;

  for ( Halfedge_const_handle lHE = aSs.halfedges_begin() ; lHE != aSs.halfedges_end() ; ++ lHE )
  {
    if ( lHE->id() > lMaxID )
      lMaxID = lHE->id() ;

    if ( !lHE->is_bisector() && handle_assigned(lHE->face()) )
      mBorders.push_back(lHE);
  }

  CGAL_POLYOFFSET_TRACE(3, "Border count: " << mBorders.size() ) ;

  CGAL_POLYOFFSET_TRACE(3, "Highest Bisector ID: " << lMaxID ) ;

  mVisitedBisectors.resize(lMaxID+1);

  ResetVisitedBisectorsMap();
}

template<class Ss, class Gt, class Cont>
typename Polygon_offset_builder_2<Ss,Gt,Cont>::Halfedge_const_handle
Polygon_offset_builder_2<Ss,Gt,Cont>::LocateHook( FT aTime, Halfedge_const_handle aBisector )
{
  CGAL_POLYOFFSET_TRACE(2,"Searching for hook at " << aTime ) ;

  Halfedge_const_handle rHook ;

  while ( aBisector->is_bisector() )
  {
    Halfedge_const_handle lPrev = aBisector->prev();
    Halfedge_const_handle lNext = aBisector->next();

    if ( !IsVisited(aBisector) )
    {
      CGAL_POLYOFFSET_TRACE(3,"Testing hook on " << e2str(*aBisector) 
                           << "\n  Next: " << e2str(*lNext) 
                           << "\n  Prev: " << e2str(*lPrev) 
                           ) ;
      
      Comparison_result lCNext = lNext->is_bisector() ? Compare_offset_against_event_time(aTime,aBisector,lNext)
                                                      : SMALLER ;

      Comparison_result lCPrev = lPrev->is_bisector() ? Compare_offset_against_event_time(aTime,lPrev,aBisector)
                                                      : SMALLER ;

      CGAL_POLYOFFSET_TRACE(3,"CPrev: " << lCPrev << " CNext: " << lCNext ) ;

      if ( lCNext != lCPrev )
      {
        CGAL_POLYOFFSET_TRACE(2, "Hook found on B" << aBisector->id() ) ;
        rHook = aBisector ;
        break ;
      }
    }
    aBisector = lNext ;
  }

  return rHook;
}

template<class Ss, class Gt, class Cont>
typename Polygon_offset_builder_2<Ss,Gt,Cont>::Halfedge_const_handle
Polygon_offset_builder_2<Ss,Gt,Cont>::LocateSeed( FT aTime )
{
  Halfedge_const_handle rSeed ;

  for ( typename Halfedge_vector::const_iterator f = mBorders.begin()
       ; f != mBorders.end() && !handle_assigned(rSeed)
       ; ++ f
      )
  {
    CGAL_POLYOFFSET_TRACE(3,"Locating hook for face E" << e2str(**f) ) ;
    rSeed = LocateHook(aTime,(*f)->next());
  }
  CGAL_POLYOFFSET_TRACE(3,"Seed found on B" << ( handle_assigned(rSeed) ? e2str(*rSeed) : "<none>" ) ) ;
  return rSeed;
}

template<class Ss, class Gt, class Cont>
void Polygon_offset_builder_2<Ss,Gt,Cont>::AddOffsetVertex( FT aTime, Halfedge_const_handle aHook, ContainerPtr aPoly )
{
  Visit(aHook);

  boost::optional<Point_2> lP = Construct_offset_point(aTime,aHook);

  if ( !lP )
    throw std::range_error("CGAL_POLYGON_OFFSET: Overflow during construction of offset vertex" ) ; // Caught by the main loop
    
  CGAL_POLYOFFSET_TRACE(1,"Constructing offset point along B" << e2str(*aHook) ) ;

  aPoly->push_back(*lP);
}

template<class Ss, class Gt, class Cont>
template<class OutputIterator>
OutputIterator Polygon_offset_builder_2<Ss,Gt,Cont>::TraceOffsetPolygon( FT aTime, Halfedge_const_handle aSeed, OutputIterator aOut )
{
  CGAL_POLYOFFSET_TRACE(1,"Tracing new offset polygon" ) ;

  ContainerPtr lPoly( new Container() ) ;

  Halfedge_const_handle lHook = aSeed ;

  AddOffsetVertex(aTime,lHook,lPoly);

  while ( true )
  {
    lHook = LocateHook(aTime,lHook->next()) ;

    if ( handle_assigned(lHook) )
    {
      if ( lHook != aSeed )
        AddOffsetVertex(aTime,lHook,lPoly);

      Halfedge_const_handle lNextBisector = lHook->opposite();

      if ( lNextBisector != aSeed && !IsVisited(lNextBisector) )
      {
        lHook = lNextBisector;
        continue;
      }
    }
    break ;
  }

  if ( lPoly->size() >= 3 )
  {
    CGAL_POLYOFFSET_TRACE(1,"Offset polygon of " << lPoly->size() << " vertices traced." ) ;
    *aOut++ = lPoly ;
  }
  else
  {
    CGAL_POLYOFFSET_TRACE(1,"Invalid offset polygon (less than 3 vertices) traced." ) ;
  }

  return aOut ;
}

template<class Ss, class Gt, class Cont>
void Polygon_offset_builder_2<Ss,Gt,Cont>::ResetVisitedBisectorsMap()
{
  std::fill(mVisitedBisectors.begin(),mVisitedBisectors.end(),0);
}

template<class Ss, class Gt, class Cont>
template<class OutputIterator>
OutputIterator Polygon_offset_builder_2<Ss,Gt,Cont>::construct_offset_contours( FT aTime, OutputIterator aOut )
{
  CGAL_precondition( aTime > static_cast<FT>(0.0) ) ;

  ResetVisitedBisectorsMap();

  CGAL_POLYOFFSET_TRACE(1,"Constructing offset polygons for offset: " << aTime ) ;
  for ( Halfedge_const_handle lSeed = LocateSeed(aTime); handle_assigned(lSeed); lSeed = LocateSeed(aTime) )
  {
    try
    {
      aOut = TraceOffsetPolygon(aTime,lSeed,aOut);
    }
    catch( std::exception const& e )
    {
      CGAL_POLYOFFSET_TRACE(0,"EXCEPTION THROWN (" << e.what() << ") during during polygon offset construction." ) ;
    }
  }  

  return aOut ;
}

CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_OFFSET_BUILDER_2_IMPL_H //
// EOF //
