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
// file          : include/CGAL/Straight_skeleton_builder_traits_2.c
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_C
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_C 1

#include <cmath>

#include <CGAL/number_utils.h>
#include <CGAL/Quotient.h>

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_SHOW_TRAITS
#  define CGAL_SSBUILDER_SHOW_TRAITS(code) code
#else
#  define CGAL_SSBUILDER_SHOW_TRAITS(code)
#endif


CGAL_BEGIN_NAMESPACE

template<class R>
boost::optional
  <
  std::pair
    <
       typename Straight_skeleton_builder_traits_2<R>::Point_2
      ,typename Straight_skeleton_builder_traits_2<R>::FT 
    >
  >
Straight_skeleton_builder_traits_2<R>::compute_event (  Segment_2 const& aA
                                                      , Segment_2 const& aB
                                                      , Segment_2 const& aC
                                                     ) const
{
  OptionalEventData rResult ;
  
  Line_2 lLineA = aA.supporting_line();
  Line_2 lLineB = aB.supporting_line();
  Line_2 lLineC = aC.supporting_line();
  
  OptionalPoint_2 lP = ConstructEventPoint(lLineA,lLineB,lLineC);
  
  if ( !!lP )
  {
    if ( IsInsideOffsetRegion(*lP,lLineA,lLineB,lLineC) )
    {
      CGAL_SSBUILDER_SHOW_TRAITS ( SS_IO_AUX::ScopedPointDrawing lDraw(*lP,3,"PotentialEvent"); )
      
      typename Rep::Compute_squared_distance_2 squared_distance;
      
      FT lTime1 = squared_distance(*lP,lLineA);
      FT lTime2 = squared_distance(*lP,lLineB);
      FT lTime3 = squared_distance(*lP,lLineC);
      FT lTime  = (lTime1 + lTime2 + lTime3 ) / 3.0 ;
      
      rResult = std::make_pair(*lP,lTime);
    }
  }
  
  return rResult ;  
}

template<class R>
Comparison_result Straight_skeleton_builder_traits_2<R>
  ::compare_events (  Segment_2 const& aAa, Segment_2 const& aAb, Segment_2 const& aAc
                    , Segment_2 const& aBa, Segment_2 const& aBb, Segment_2 const& aBc 
                   ) const 
{
  OptionalEventData lEventA = compute_event(aAa,aAb,aAc);
  OptionalEventData lEventB = compute_event(aBa,aBb,aBc);
  CGAL_assertion( !!lEventA ) ;
  CGAL_assertion( !!lEventB ) ;
  FT lTimeA = lEventA->second ;
  FT lTimeB = lEventB->second;
  return compare(lTimeA,lTimeB);
}                   

template<class R>
bool Straight_skeleton_builder_traits_2<R>::is_event_inside_bounded_offset_zone 
                                            (  Segment_2 const& aA
                                             , Segment_2 const& aB
                                             , Segment_2 const& aC
                                             , Segment_2 const& aEdge
                                             , Segment_2 const& aEdgeLeft
                                             , Segment_2 const& aEdgeRight
                                            ) const  
{
  OptionalEventData lEvent = compute_event(aA,aB,aC);
  CGAL_assertion( !!lEvent ) ;
  
  Line_2 lEdgeLeftLine  = aEdgeLeft .supporting_line();
  Line_2 lEdgeLine      = aEdge     .supporting_line();
  Line_2 lEdgeRightLine = aEdgeRight.supporting_line();
  
  Line_2 lBisectorL = ConstructBisector(lEdgeLeftLine,lEdgeLine);
  Line_2 lBisectorR = ConstructBisector(lEdgeLine,lEdgeRightLine);
  
  Point_2 lP = lEvent->first;
  
  typename Rep::Oriented_side_2 oriented_side ;
  
  return    oriented_side(lBisectorL,lP) == CGAL::ON_NEGATIVE_SIDE 
         && oriented_side(lEdgeLine ,lP) == CGAL::ON_POSITIVE_SIDE 
         && oriented_side(lBisectorR,lP) == CGAL::ON_POSITIVE_SIDE ;
}                             
                             

template<class R>
typename Straight_skeleton_builder_traits_2<R>::OptionalPoint_2
Straight_skeleton_builder_traits_2<R>::Intersect ( Line_2 const& aA, Line_2 const& aB ) const
{
  typename Rep::Compute_squared_distance_2 squared_distance ;
  
  OptionalPoint_2 rResult ; 
  Object lIsec = Rep().intersect_2_object()(aA,aB);
  if ( !lIsec.is_empty() )
  {
    Point_2 lP ; 
    if ( assign(lP,lIsec) )
    {
      //
      // Hack for unfiltered floating-point kernels!!!!!!!
      //
      FT const lDiskR = 1e-5 ;
      FT lDevL = squared_distance(aA,lP);
      FT lDevR = squared_distance(aB,lP);
      if ( lDevL <= lDiskR && lDevR <= lDiskR )
        rResult.reset(lP);
    }  
  }
  return rResult ;
}

template<class R>
typename Straight_skeleton_builder_traits_2<R>::Line_2
Straight_skeleton_builder_traits_2<R>::ConstructBisector(  Line_2 const& aA
                                                         , Line_2 const& aB
                                                        ) const
{
  Point_2  lBaseP ;
  Vector_2 lDir ;
  
  OptionalPoint_2 lI = Intersect(aA,aB);
  if ( lI )
  {
    lBaseP = *lI ;
    
    Direction_2 lDirA = -aA.direction() ;
    Direction_2 lDirB =  aB.direction();
    
    double lAngleA = std::atan2(to_double(lDirA.dy()), to_double(lDirA.dx()));
    double lAngleB = std::atan2(to_double(lDirB.dy()), to_double(lDirB.dx()));
    
    double const cPi     = CGAL_PI ;
    double const cTwoPi  = CGAL_PI * 2 ;
    double const cFourPi = CGAL_PI * 4 ;
    
    double lNormAngleA = lAngleA < 0 ? lAngleA + cTwoPi : lAngleA ;
    double lNormAngleB = lAngleB < 0 ? lAngleB + cTwoPi : lAngleB ;
     
    double lSweep = lNormAngleA == lNormAngleB 
                      ? 0 : fmod((lNormAngleA-lNormAngleB) + cFourPi,cTwoPi) ;
                               
    double lPhi = lSweep / 2.0 ;

    FT s = std::sin(lPhi) ;
    FT c = std::cos(lPhi) ;
  
    typename Rep::Aff_transformation_2 Rot( CGAL::Rotation(), s ,c ) ;
    
    lDir = lDirB.transform(Rot).to_vector();
  }
  else
  {
    Point_2 lQA = aA.point(0);
    
    typename Rep::Compute_squared_distance_2 squared_distance ;
    typename Rep::Compute_squared_length_2   squared_length ;
    
    FT lDist = squared_distance(aA,aB);
    
    lDir = aA.to_vector();
    
    Vector_2 lN = lDir.perpendicular( static_cast<Orientation>(aA.oriented_side(aB.point(0))) ) ;
                                    
    FT lLen = CGAL::sqrt( squared_length(lN) ) ;
    
    lBaseP = lQA + Vector_2(lN.x()*lDist,lN.y()*lDist,lLen);
  }
  
  return Line_2(lBaseP,lDir);
}

template<class R>
typename Straight_skeleton_builder_traits_2<R>::OptionalPoint_2
Straight_skeleton_builder_traits_2<R>::ConstructEventPoint (  Line_2 const& aA
                                                            , Line_2 const& aB
                                                            , Line_2 const& aC
                                                           ) const
{
  Line_2 const& lBisecL = ConstructBisector(aA,aB);
  Line_2 const& lBisecR = ConstructBisector(aB,aC);
  
  /*
  CGAL_SSBUILDER_SHOW_TRAITS
  (
    SS_IO_AUX::ScopedSegmentDrawing lDrawL(lBisecL,0,"UnboundBisector");
    SS_IO_AUX::ScopedSegmentDrawing lDrawR(lBisecR,0,"UnboundBisector");
  )
  */
  
  return Intersect(lBisecL,lBisecR);
}                                                            

template<class R>
bool Straight_skeleton_builder_traits_2<R>::IsInsideOffsetRegion (  Point_2 const& aP
                                                                  , Line_2  const& aLineA
                                                                  , Line_2  const& aLineB
                                                                  , Line_2  const& aLineC
                                                                 ) const 
{
  typename Rep::Oriented_side_2 oriented_side ;
  
  return    oriented_side(aLineA,aP) == CGAL::ON_POSITIVE_SIDE 
         && oriented_side(aLineB,aP) == CGAL::ON_POSITIVE_SIDE 
         && oriented_side(aLineC,aP) == CGAL::ON_POSITIVE_SIDE ;
}

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_C //
// EOF //
