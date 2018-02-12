// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_POLYGON_OFFSET_BUILDER_2_H
#define CGAL_POLYGON_OFFSET_BUILDER_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/disable_warnings.h>

#include <vector>
#include <algorithm>

#include <boost/shared_ptr.hpp>
#include <boost/optional/optional.hpp>

#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>

#include <CGAL/Polygon_offset_builder_traits_2.h>

namespace CGAL {

template<class Traits_, class SSkel_>
struct Default_polygon_offset_builder_2_visitor
{
  typedef Traits_ Traits ;
  typedef SSkel_  SSkel ;

  typedef typename SSkel::Halfedge_const_handle Halfedge_const_handle ;
  typedef typename SSkel::Vertex_const_handle   Vertex_const_handle ;
  
  typedef typename Traits::FT      FT ;
  typedef typename Traits::Point_2 Point_2 ;
  
  void on_construction_started ( FT ) const {}
  
  void on_offset_contour_started() const {}
  
  void on_offset_point ( Point_2 const& ) const {}

  Point_2 on_offset_point_overflowed( Halfedge_const_handle aBisector ) const
  {
    ::CGAL::warning_fail( "", __FILE__, __LINE__, "! Unable to compute polygon offset point due to overflow !" ) ;
    return aBisector->vertex()->point() ;
  }
  
  void on_offset_contour_finished ( bool /*is_complete*/ ) const {}
  
  void on_construction_finished () const {}
  
  void on_error( char const* ) const {}
} ;

template<class Ss_, class Traits_, class Container_, class Visitor_ = Default_polygon_offset_builder_2_visitor<Traits_,Ss_> >
class Polygon_offset_builder_2
{
public :

  typedef Ss_        Ss ;
  typedef Traits_    Traits ;
  typedef Container_ Container ;
  typedef Visitor_   Visitor ;

  typedef typename Traits::Kernel K ;
  
  typedef typename Traits::FT      FT ;
  typedef typename Traits::Point_2 Point_2 ;
  
  typedef boost::optional<Point_2> OptionalPoint_2 ;

  typedef boost::shared_ptr<Container> ContainerPtr ;

  Polygon_offset_builder_2( Ss const& aSs, Traits const& aTraits = Traits(), Visitor const& aVisitor = Visitor() )  ;

  template<class OutputIterator>
  OutputIterator construct_offset_contours( FT aTime, OutputIterator aOut ) ;

  struct Bisector_data
  {
    Bisector_data() : IsVisited(false), IsUsedSeed(false) {}
    
    bool IsVisited ;
    bool IsUsedSeed ;
  } ;
  
private:

  typedef typename Ss::Vertex                Vertex ;
  typedef typename Ss::Halfedge_handle       Halfedge_handle  ;
  typedef typename Ss::Halfedge_const_handle Halfedge_const_handle  ;
  typedef typename Ss::Vertex_const_handle   Vertex_const_handle  ;
  
  typedef std::vector<Halfedge_const_handle> Halfedge_vector ;

  typedef typename Traits::Segment_2        Segment_2 ;
  typedef typename Traits::Trisegment_2     Trisegment_2 ;
  typedef typename Traits::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef CGAL_SS_i::Triedge<Halfedge_handle> Triedge ;

  enum Hook_position { SOURCE, TARGET, INSIDE } ;
  
  static char const* Hook_position2Str( Hook_position aPos )
  {
    return aPos == SOURCE ? "source" : ( aPos == TARGET ? "target" : "inside" ) ;
  }
  
  Halfedge_const_handle LocateHook( FT aTime, Halfedge_const_handle aBisector, bool aIncludeLastBisector, Hook_position& rPos ) ;

  void AddOffsetVertex( FT aTime, Halfedge_const_handle aHook, ContainerPtr aPoly ) ;

  template<class OutputIterator>
  OutputIterator TraceOffsetPolygon( FT aTime, Halfedge_const_handle aHook, OutputIterator aOut ) ;

  Halfedge_const_handle LocateSeed( FT aTime, Halfedge_const_handle aBorder ) ;
  
  Halfedge_const_handle LocateSeed( FT aTime ) ;

  Bisector_data const& GetBisectorData ( Halfedge_const_handle aBisector ) const { return mBisectorData[aBisector->id()] ; }
  Bisector_data&       GetBisectorData ( Halfedge_const_handle aBisector )       { return mBisectorData[aBisector->id()] ; }
    
  bool IsVisited( Halfedge_const_handle aBisector ) { return GetBisectorData(aBisector).IsVisited ; }

  bool IsUsedSeed( Halfedge_const_handle aBisector ) { return GetBisectorData(aBisector).IsUsedSeed ; }
  
  void Visit( Halfedge_const_handle aBisector ) { GetBisectorData(aBisector).IsVisited = true ; }
  
  void SetIsUsedSeed( Halfedge_const_handle aBisector ) { GetBisectorData(aBisector).IsUsedSeed = true ; }

  inline Segment_2 CreateSegment ( Halfedge_const_handle aH ) const
  {
    Point_2 s = aH->opposite()->vertex()->point() ;
    Point_2 t = aH->vertex()->point() ;
    return K().construct_segment_2_object()(s,t);
  }

  Trisegment_2_ptr CreateTrisegment ( Triedge const& aTriedge ) const
  {
    CGAL_precondition( aTriedge.is_valid() ) ;
    
    if ( aTriedge.is_skeleton() )
    {
      return Construct_ss_trisegment_2(mTraits)(CreateSegment(aTriedge.e0())
                                               ,CreateSegment(aTriedge.e1())
                                               ,CreateSegment(aTriedge.e2())
                                               );
    }
    else 
    {
      return Trisegment_2_ptr() ;
    }
  }
  
  Trisegment_2_ptr CreateTrisegment ( Vertex_const_handle aNode ) const ;
  
  Vertex_const_handle GetSeedVertex ( Vertex_const_handle   aNode
                                    , Halfedge_const_handle aBisector
                                    , Halfedge_const_handle aEa
                                    , Halfedge_const_handle aEb 
                                    ) const ;

  bool Is_bisector_defined_by ( Halfedge_const_handle aBisector, Halfedge_const_handle aEa, Halfedge_const_handle aEb ) const 
  { 
    return    ( aBisector->defining_contour_edge() == aEa && aBisector->opposite()->defining_contour_edge() == aEb ) 
           || ( aBisector->defining_contour_edge() == aEb && aBisector->opposite()->defining_contour_edge() == aEa ) ;  
  }
  
  Comparison_result Compare_offset_against_event_time( FT aT, Vertex_const_handle aNode ) const
  {
    CGAL_precondition( aNode->is_skeleton() ) ;
    
    Comparison_result r = aNode->has_infinite_time() ? SMALLER
                                                     : static_cast<Comparison_result>(Compare_offset_against_event_time_2(mTraits)(aT,CreateTrisegment(aNode)));
    
    return r ;
  }
public:
  boost::optional<Point_2> Construct_offset_point( FT aT, Halfedge_const_handle aBisector ) const
  {
    CGAL_assertion(aBisector->is_bisector());
    CGAL_assertion(handle_assigned(aBisector->opposite()));
    CGAL_assertion(aBisector->opposite()->is_bisector());

    Halfedge_const_handle lBorderA = aBisector->defining_contour_edge();
    Halfedge_const_handle lBorderB = aBisector->opposite()->defining_contour_edge();
    
    Vertex_const_handle lNodeS = aBisector->opposite()->vertex();
    Vertex_const_handle lNodeT = aBisector->vertex();
    
    
    // If aBisector is not a border bisector the offset point construction needs to get to seed event
    Trisegment_2_ptr lSeedEvent ;
    if ( aBisector->is_inner_bisector() )
    {
      CGAL_assertion ( lNodeS->is_skeleton() ) ;
      CGAL_assertion ( lNodeT->is_skeleton() ) ;
      
      Vertex_const_handle lSeedNode = aBisector->slope() == POSITIVE ? lNodeS : lNodeT ;
      
      lSeedEvent = CreateTrisegment(lSeedNode) ;
      
      CGAL_POLYOFFSET_TRACE(3,"Seed node for " << e2str(*aBisector) << " is " << v2str(*lSeedNode) << " event=" << lSeedEvent ) ;
    }

    OptionalPoint_2 p = Construct_offset_point_2(mTraits)(aT
                                                         ,CreateSegment(lBorderA)
                                                         ,CreateSegment(lBorderB)
                                                         ,lSeedEvent
                                                         );
    CGAL_stskel_intrinsic_test_assertion
    ( 
      !p 
      || 
      ( p && !CGAL_SS_i::is_possibly_inexact_distance_clearly_not_zero
              ( CGAL_SS_i::squared_distance_from_point_to_lineC2(p->x()
                                                                ,p->y()
                                                                ,lNodeS->point().x()
                                                                ,lNodeS->point().y()
                                                                ,lNodeT->point().x()
                                                                ,lNodeT->point().y()
                                                                ).to_nt()
              )
      )
    ) ;
    
    return p ;
  }
private:
  void ResetBisectorData();

  Traits const&              mTraits ;
  Visitor const&             mVisitor ;
  Halfedge_vector            mBorders ;
  std::vector<Bisector_data> mBisectorData;
  OptionalPoint_2            mLastPoint ; 
  CGAL_POLYOFFSET_DEBUG_CODE( int mStepID ; )
};

} // end namespace CGAL

#include <CGAL/Straight_skeleton_2/Polygon_offset_builder_2_impl.h>

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_OFFSET_BUILDER_2_H //
// EOF //



