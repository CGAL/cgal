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
// file          : include/CGAL/Straight_skeleton_builder_traits_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H 1

#include <algorithm>
#include <boost/optional.hpp>

#include <CGAL/Straight_skeleton_aux.h>

CGAL_BEGIN_NAMESPACE

template<class R>
class Straight_skeleton_builder_traits_2 
{
public:
  
  typedef R Rep ;

  typedef typename Rep::Point_2 Point_2 ;
  typedef typename Rep::FT      FT ;
  
public:

  typedef std::pair<Point_2,Point_2> Point_2_Pair ;
  
  typedef std::pair<Point_2,FT> EventData ;
  
  typedef boost::optional<EventData> OptionalEventData ;
  
  OptionalEventData compute_event (  Point_2_Pair const& aA
                                   , Point_2_Pair const& aB
                                   , Point_2_Pair const& aC
                                  ) const ;
  
  Comparison_result compare_events (  Point_2_Pair const& aLA
                                    , Point_2_Pair const& aLB
                                    , Point_2_Pair const& aLC
                                    , Point_2_Pair const& aRA
                                    , Point_2_Pair const& aRB
                                    , Point_2_Pair const& aRC
                                   ) const ;
                                   
  Comparison_result compare_events_distance_to_seed ( Point_2             aSeedP
                                                    , Point_2_Pair const& aLA
                                                    , Point_2_Pair const& aLB
                                                    , Point_2_Pair const& aLC
                                                    , Point_2_Pair const& aRA
                                                    , Point_2_Pair const& aRB
                                                    , Point_2_Pair const& aRC
                                                  ) const ;
                                   
  Comparison_result compare_events_distance_to_seed ( Point_2_Pair const& aSeedA
                                                    , Point_2_Pair const& aSeedB
                                                    , Point_2_Pair const& aSeedC
                                                    , Point_2_Pair const& aLA
                                                    , Point_2_Pair const& aLB
                                                    , Point_2_Pair const& aLC
                                                    , Point_2_Pair const& aRA
                                                    , Point_2_Pair const& aRB
                                                    , Point_2_Pair const& aRC
                                                  ) const ;
                                                  
  bool is_event_inside_bounded_offset_zone (  Point_2_Pair const& aA
                                            , Point_2_Pair const& aB
                                            , Point_2_Pair const& aC
                                            , Point_2_Pair const& aEdge
                                            , Point_2_Pair const& aEdgeLeft
                                            , Point_2_Pair const& aEdgeRight
                                           ) const ; 

private :

  typedef typename Rep::Vector_2    Vector_2 ;
  typedef typename Rep::Direction_2 Direction_2 ;
  typedef typename Rep::Line_2      Line_2 ;
  typedef typename Rep::Segment_2   Segment_2 ;
  
  typename Rep::Construct_line_2 construct_line ;
  
  typedef boost::optional<Point_2> OptionalPoint_2 ;
  
  OptionalPoint_2 Intersect ( Line_2 const& aA, Line_2 const& aB ) const ;
  
  Line_2 ConstructBisector( Line_2 const& aA, Line_2 const& aB ) const ;
  
  OptionalPoint_2 ConstructEventPoint ( Line_2 const& aA, Line_2 const& aB, Line_2 const& aC ) const ;
 
  bool IsInsideOffsetRegion (  Point_2   const& aP
                             , Line_2 const& aLineA
                             , Line_2 const& aLineB
                             , Line_2 const& aLineC
                           ) const ;
                           
} ;

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Straight_skeleton_builder_traits_2.C>
#endif


#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H //
// EOF //
