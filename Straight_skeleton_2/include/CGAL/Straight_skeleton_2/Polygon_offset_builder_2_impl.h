// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_POLYGON_OFFSET_BUILDER_2_IMPL_H
#define CGAL_POLYGON_OFFSET_BUILDER_2_IMPL_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <algorithm>
#include <vector>

namespace CGAL {


template<class Ss, class Gt, class Cont, class Visitor>
Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::Polygon_offset_builder_2( Ss const& aSs, Traits const& aTraits, Visitor const& aVisitor )
  :
   mTraits (aTraits)
  ,mVisitor(aVisitor)
{

  int lMaxID = -1 ;

  for ( Halfedge_const_handle lHE = aSs.halfedges_begin() ; lHE != aSs.halfedges_end() ; ++ lHE )
  {
    if ( lHE->id() > lMaxID )
      lMaxID = lHE->id() ;

    if ( !lHE->is_bisector() && handle_assigned(lHE->face()) )
      mBorders.push_back(lHE);
  }

  CGAL_POLYOFFSET_TRACE(2, "Border count: " << mBorders.size() ) ;

  CGAL_POLYOFFSET_TRACE(2, "Highest Bisector ID: " << lMaxID ) ;

  mBisectorData.resize(lMaxID+1);

  ResetBisectorData();
}


template<class Ss, class Gt, class Cont, class Visitor>
typename Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::Halfedge_const_handle
Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::LocateHook( FT                    aTime
                                                        , Halfedge_const_handle aBisector
                                                        , bool                  aIncludeLastBisector
                                                        , Hook_position&        rPos
                                                        )
{
  CGAL_POLYOFFSET_TRACE(2,"Locate hook at " << aTime ) ;
  CGAL_POLYOFFSET_TRACE(2,"Start halfedge: " << e2str(*aBisector) ) ;

  Halfedge_const_handle rHook ;

  while ( aBisector->is_bisector() && ( aIncludeLastBisector ? true : aBisector->prev()->is_bisector() ) )
  {
    Halfedge_const_handle lPrev = aBisector->prev();
    Halfedge_const_handle lNext = aBisector->next();

    CGAL_POLYOFFSET_TRACE(2,"Testing hook on " << e2str(*aBisector) ) ;
    CGAL_POLYOFFSET_TRACE(4, "Next: " << e2str(*lNext) << " - Prev: " << e2str(*lPrev) ) ;

    if ( !IsVisited(aBisector) )
    {
      if ( aBisector->slope() != ZERO )
      {
        // A hook is found here if 'aTime' is within the bisector time interval.
        //
        // Depending on the bisector slope, src-time might be smaller or larger than tgt-time,
        // so the test is:
        //
        //  (src-time <= time <= tgt-time ) || ( tgt-time <= time <= src-time )
        //

        Comparison_result lTimeWrtSrcTime = lPrev->is_bisector() ? Compare_offset_against_event_time(aTime,lPrev    ->vertex()) : LARGER ;
        Comparison_result lTimeWrtTgtTime = lNext->is_bisector() ? Compare_offset_against_event_time(aTime,aBisector->vertex()) : LARGER ;
        CGAL_POLYOFFSET_TRACE(3,"  TimeWrtSrcTime: " << lTimeWrtSrcTime << " TimeWrtTgtTime: " << lTimeWrtTgtTime ) ;

        //
        // The above test expressed in terms of comparisons of src/tgt time against aTime is:
        //
        // (     ( time-wrt-src-time == ZERO || time-wrt-src-time == SMALLER )
        //    && ( time-wrt-tgt-time == ZERO || time-wrt-tgt-time == LARGER  )
        // )
        //
        // || -the same with src/tgt inverted-
        //
        // (      ( time-wrt-tgt-time == ZERO || time-wrt-tgt-time == SMALLER )
        //     && ( time-wrt-src-time == ZERO || time-wrt-src-time == LARGER  )
        // )
        //
        // But since bisectors of slope zero are skipped, both comparisons cannot be zero, thus, the test above is really:
        //
        //    ( ( time-wrt-src-time == ZERO || time-wrt-src-time == SMALLER )  && (                              time-wrt-tgt-time == LARGER ) )
        // || ( (                              time-wrt-src-time == SMALLER )  && ( time-wrt-tgt-time == ZERO || time-wrt-tgt-time == LARGER ) )
        // || ( ( time-wrt-tgt-time == ZERO || time-wrt-tgt-time == SMALLER )  && (                              time-wrt-src-time == LARGER ) )
        // || ( (                              time-wrt-tgt-time == SMALLER )  && ( time-wrt-src-time == ZERO || time-wrt-src-time == LARGER ) )
        //
        // Which actually boils down to this:
        //
        if ( lTimeWrtSrcTime != lTimeWrtTgtTime )
        {
          CGAL_stskel_intrinsic_test_assertion( !CGAL_SS_i::is_time_clearly_not_within_possibly_inexact_bisector_time_interval(aTime,aBisector) ) ;

          rPos = ( lTimeWrtTgtTime == EQUAL ? TARGET : lTimeWrtSrcTime == EQUAL ? SOURCE : INSIDE ) ;

          rHook = aBisector ;

          CGAL_POLYOFFSET_TRACE(2, "  Hook found here at " << Hook_position2Str(rPos) ) ;

          break ;
        }
        else
        {
          CGAL_stskel_intrinsic_test_assertion( !CGAL_SS_i::is_time_clearly_within_possibly_inexact_bisector_time_interval(aTime,aBisector) ) ;

          CGAL_POLYOFFSET_TRACE(2, "  Hook not found here.") ;
        }
      }
      else
      {
        CGAL_POLYOFFSET_TRACE(2,"Bisector is a roof peak (zero slope).");
      }
    }
    else
    {
      CGAL_POLYOFFSET_TRACE(2,"Bisector already visited");
    }
    aBisector = lPrev ;
  }

  return rHook;
}

template<class Ss, class Gt, class Cont, class Visitor>
typename Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::Halfedge_const_handle
Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::LocateSeed( FT aTime, Halfedge_const_handle aBorder )
{
  CGAL_POLYOFFSET_TRACE(2,"\nSearching for a starting seed in face " << e2str(*aBorder) ) ;

  Hook_position lPos ;
  Halfedge_const_handle rSeed = LocateHook(aTime,aBorder->prev(),false,lPos);
  if ( handle_assigned(rSeed) )
  {
    if ( !IsUsedSeed(rSeed) )
    {
      SetIsUsedSeed(rSeed);

      CGAL_postcondition( handle_assigned(rSeed->prev()) && rSeed->prev()->is_bisector() ) ;

      // If a seed hook is found at a bisector's source,
      // the next hook will be found at the previous bisector's target, which would be a mistake.
      // So, we modify the seed to be the target() of the previous halfedge instead.
      if ( lPos == SOURCE )
        rSeed = rSeed->prev() ;
      CGAL_POLYOFFSET_TRACE(2,"Pos at source switched to pos at target on " << e2str(*rSeed) ) ;
    }
    else
    {
      CGAL_POLYOFFSET_TRACE(2,"Seed already used. Discarded");
      rSeed = Halfedge_const_handle();
    }
  }
  return rSeed ;
}


template<class Ss, class Gt, class Cont, class Visitor>
typename Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::Halfedge_const_handle
Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::LocateSeed( FT aTime )
{
  CGAL_POLYOFFSET_TRACE(2,"Searching for a starting seed at " << aTime ) ;

  Halfedge_const_handle rSeed ;

  for ( typename Halfedge_vector::const_iterator f = mBorders.begin()
       ; f != mBorders.end() && !handle_assigned(rSeed)
       ; ++ f
      )
    rSeed = LocateSeed(aTime,*f);

  CGAL_POLYOFFSET_TRACE(2,"Found seed: " << eh2str(rSeed) ) ;

  return rSeed;
}


template<class Ss, class Gt, class Cont, class Visitor>
void Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::AddOffsetVertex( FT                    aTime
                                                                  , Halfedge_const_handle aHook
                                                                  , ContainerPtr          aPoly
                                                                  )
{
  OptionalPoint_2 lP = Construct_offset_point(aTime,aHook);

  if ( !lP )
    lP = mVisitor.on_offset_point_overflowed(aHook) ;

  CGAL_postcondition(bool(lP));

  CGAL_POLYOFFSET_TRACE(1,"Found offset point p=" << p2str(*lP) << " at offset " << aTime << " along bisector " << e2str(*aHook) << " reaching " << v2str(*aHook->vertex()) ) ;

  mVisitor.on_offset_point(*lP);

  if ( lP != mLastPoint )
  {
    aPoly->push_back(*lP);
    mLastPoint = lP ;
  }
  else
  {
    CGAL_POLYOFFSET_TRACE(1,"Duplicate point. Ignored");
  }

  CGAL_POLYOFFSET_DEBUG_CODE( ++ mStepID ) ;
}

template<class Ss, class Gt, class Cont, class Visitor>
template<class OutputIterator>
OutputIterator Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::TraceOffsetPolygon( FT aTime, Halfedge_const_handle aSeed, OutputIterator aOut )
{
  CGAL_POLYOFFSET_TRACE(1,"\nTracing new offset polygon" ) ;

  ContainerPtr lPoly( new Container() ) ;

  mVisitor.on_offset_contour_started();

  Halfedge_const_handle lHook = aSeed ;
  std::vector <Halfedge_const_handle > visited_hooks;
  do
  {
    CGAL_POLYOFFSET_TRACE(1,"STEP " << mStepID ) ;

    Halfedge_const_handle lLastHook = lHook ;

    Hook_position lPos ;
    lHook = LocateHook(aTime,lHook->prev(),true,lPos) ;

    if ( handle_assigned(lHook) )
    {
      CGAL_POLYOFFSET_TRACE(4, "returned Hook: " << e2str(*lHook));

      CGAL_assertion( lHook->slope() != ZERO );
      AddOffsetVertex(aTime,lHook, lPoly);

      Visit(lHook);
      visited_hooks.push_back(lHook);

      CGAL_POLYOFFSET_TRACE(2,"Marking hook, B" << lHook->id() << ", as visited." ) ;

      lHook = lHook->opposite();
    }

    Visit(lLastHook);
    visited_hooks.push_back(lLastHook);

    CGAL_POLYOFFSET_TRACE(2,"Marking last hook, B" << lLastHook->id() << ", as visited." ) ;
  }
  while ( handle_assigned(lHook) && lHook != aSeed && !IsVisited(lHook)) ;

  bool lComplete = ( lHook == aSeed )  ;

  CGAL_POLYOFFSET_TRACE(1,"Offset polygon of " << lPoly->size() << " vertices traced." << ( lComplete ? "COMPLETE" : "INCOMPLETE" ) ) ;

  // On paper, lComplete == true should imply that lPoly->size() >= 3, but since the constructions
  // might not be exact, you can have cases where the offset points are actually duplicates
  // and so the end polygon has size < 3. It is ignored in that case.
  if ( lComplete && lPoly->size() < 3 )
    lComplete = false;

  mVisitor.on_offset_contour_finished( lComplete );

  if ( lComplete )
  {
    *aOut++ = lPoly ;
  }
  else
  {
    for (std::size_t k=0;k<visited_hooks.size();++k)
    {
      GetBisectorData( visited_hooks[k] ).IsVisited=false;
    }
  }

  return aOut ;
}

template<class Ss, class Gt, class Cont, class Visitor>
void Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::ResetBisectorData()
{
  std::fill(mBisectorData.begin(),mBisectorData.end(), Bisector_data() );
}

template<class Ss, class Gt, class Cont, class Visitor>
template<class OutputIterator>
OutputIterator Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::construct_offset_contours( FT aTime, OutputIterator aOut )
{
  CGAL_precondition( aTime > static_cast<FT>(0) ) ;

  CGAL_POLYOFFSET_DEBUG_CODE( mStepID = 0 ) ;

  mVisitor.on_construction_started(aTime);

  mLastPoint = boost::none ;

  ResetBisectorData();

  CGAL_POLYOFFSET_TRACE(1,"Constructing offset polygons for offset: " << aTime ) ;
  for ( Halfedge_const_handle lSeed = LocateSeed(aTime); handle_assigned(lSeed); lSeed = LocateSeed(aTime) )
    aOut = TraceOffsetPolygon(aTime,lSeed,aOut);

  mVisitor.on_construction_finished();

  return aOut ;
}

template<class Ss, class Gt, class Cont, class Visitor>
typename Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::Trisegment_2_ptr
Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::CreateTrisegment ( Vertex_const_handle aNode ) const
{
  CGAL_precondition(handle_assigned(aNode));

  Trisegment_2_ptr r ;

  CGAL_POLYOFFSET_TRACE(3,"Creating Trisegment for " << v2str(*aNode) ) ;

  if ( aNode->is_skeleton() )
  {
    Triedge const& lEventTriedge =  aNode->event_triedge() ;

    r = CreateTrisegment(lEventTriedge) ;

    CGAL_stskel_intrinsic_test_assertion
    (
      !CGAL_SS_i::is_possibly_inexact_distance_clearly_not_equal_to( Construct_ss_event_time_and_point_2(mTraits)(r)->get<0>()
                                                                   , aNode->time()
                                                                   )
    ) ;

    CGAL_POLYOFFSET_TRACE(3,"Event triedge=" << lEventTriedge ) ;

    if ( r->degenerate_seed_id() == Trisegment_2::LEFT )
    {
     CGAL_POLYOFFSET_TRACE(3,"Left seed is degenerate." ) ;

      Vertex_const_handle lLeftSeed = GetSeedVertex(aNode
                                                   ,aNode->primary_bisector()->prev()->opposite()
                                                   ,lEventTriedge.e0()
                                                   ,lEventTriedge.e1()
                                                   ) ;
      if ( handle_assigned(lLeftSeed) )
        r->set_child_l( CreateTrisegment(lLeftSeed) ) ; // Recursive call
    }
    else if ( ! aNode->is_split() && r->degenerate_seed_id() == Trisegment_2::RIGHT )
    {
      CGAL_POLYOFFSET_TRACE(3,"Right seed is degenerate." ) ;

      Vertex_const_handle lRightSeed = GetSeedVertex(aNode
                                                    ,aNode->primary_bisector()->opposite()->next()
                                                    ,lEventTriedge.e1()
                                                    ,lEventTriedge.e2()
                                                    ) ;
      if ( handle_assigned(lRightSeed) )
        r->set_child_r( CreateTrisegment(lRightSeed) ) ; // Recursive call
    }
  }

  return r ;
}

template<class Ss, class Gt, class Cont, class Visitor>
typename Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::Vertex_const_handle
Polygon_offset_builder_2<Ss,Gt,Cont,Visitor>::GetSeedVertex ( Vertex_const_handle   aNode
                                                            , Halfedge_const_handle aBisector
                                                            , Halfedge_const_handle aEa
                                                            , Halfedge_const_handle aEb
                                                            ) const
{
  Vertex_const_handle rSeed ;

  if ( Is_bisector_defined_by(aBisector,aEa,aEb) )
  {
    rSeed = aBisector->vertex();

    CGAL_POLYOFFSET_TRACE(3,"Seed of N" << aNode->id() << " for vertex (E" << aEa->id() << ",E" << aEb->id() << ") directly found: " << v2str(*rSeed) ) ;
  }
  else
  {
    typedef typename Vertex::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator ;

    Halfedge_around_vertex_const_circulator cb = aNode->halfedge_around_vertex_begin() ;
    Halfedge_around_vertex_const_circulator c  = cb ;
    do
    {
      Halfedge_const_handle lBisector = *c ;
      if ( Is_bisector_defined_by(lBisector,aEa,aEb) )
      {
        rSeed = lBisector->opposite()->vertex();
        CGAL_POLYOFFSET_TRACE(3,"Seed of N" << aNode->id() << " for vertex (E" << aEa->id() << ",E" << aEb->id() << ") indirectly found: V" << rSeed->id() ) ;
      }
    }
    while ( !handle_assigned(rSeed) && ++ c != cb ) ;
  }

  return rSeed ;
}

} // end namespace CGAL

#endif // CGAL_POLYGON_OFFSET_BUILDER_2_IMPL_H //
// EOF //
