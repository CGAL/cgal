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
// file          : include/CGAL/Polygon_offset_builder_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================

#ifndef CGAL_POLYGON_OFFSET_BUILDER_2_H
#define CGAL_POLYGON_OFFSET_BUILDER_2_H 1

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>

#include <CGAL/Polygon_offset_builder_traits_2.h>

CGAL_BEGIN_NAMESPACE

template<class Sls_, class Traits_, class Polygon_2_>
class Polygon_offset_builder_2
{
public :

  typedef Sls_       Sls ;
  typedef Traits_    Traits ;
  typedef Polygon_2_ Polygon_2 ;

  typedef typename Traits::FT FT ;

  typedef typename Polygon_2::Point_2 Point_2 ;

  typedef boost::shared_ptr<Polygon_2> Polygon_2_Ptr ;

  Polygon_offset_builder_2( Sls const& aSls, Traits const& aTraits = Traits() )  ;

  template<class OutputIterator>
  OutputIterator Create ( FT aTime, OutputIterator aOut ) ;

private:

  typedef typename Sls::Halfedge_const_handle Halfedge_const_handle  ;

  typedef std::vector<Halfedge_const_handle> Halfedge_vector ;

  typedef tuple<Point_2,Point_2> Edge ;

  typedef tuple<Edge,Edge,Edge> Edge_triple ;

  bool handled_assigned( Halfedge_const_handle aH ) const
  {
    const Halfedge_const_handle cNull ;
    return aH != cNull ;
  }

  Halfedge_const_handle LocateHook( FT aTime, Halfedge_const_handle aBisector, bool aAbove ) ;

  template<class OutputIterator>
  OutputIterator TraceOffsetPolygon( FT aTime, Halfedge_const_handle aHook, OutputIterator aOut ) ;

  Halfedge_const_handle LocateSeed( FT aTime ) ;

  bool IsVisited( Halfedge_const_handle aBisector ) { return mVisitedBisectors[aBisector->id()]; }

  void Visit( Halfedge_const_handle aBisector ) { mVisitedBisectors[aBisector->id()] = true ; }

  static inline Edge GetEdge ( Halfedge_const_handle aH )
  {
    return make_tuple(aH->opposite()->vertex()->point(),aH->vertex()->point());
  }

  static inline Edge_triple GetEdgeTriple ( Halfedge_const_handle aE0
                                          , Halfedge_const_handle aE1
                                          , Halfedge_const_handle aE2
                                          )
  {
    return make_tuple(GetEdge(aE0),GetEdge(aE1),GetEdge(aE2));
  }

  Comparison_result Compare_offset_against_event_time( FT aT, Halfedge_const_handle aBisector, Halfedge_const_handle aNextBisector ) const
  {
    CGAL_assertion(aBisector->is_bisector());
    CGAL_assertion(handle_assigned(aBisector->opposite()));
    CGAL_assertion(aBisector->opposite()->is_bisector());
    CGAL_assertion(aNextBisector->is_bisector());
    CGAL_assertion(handle_assigned(aNextBisector->opposite()));
    CGAL_assertion(aNextBisector->opposite()->is_bisector());

    Halfedge_const_handle lBorderA = aBisector->defining_border();
    Halfedge_const_handle lBorderB = aBisector->opposite()->defining_border();
    Halfedge_const_handle lBorderC = aNextBisector->opposite()->defining_border();

    return Compare_offset_against_event_time_2<Traits>(mTraits)()(aT,GetEdgeTriple(lBorderA,lBorderB,lBorderC));
  }

  Point_2 Construct_offset_point( FT aT, Halfedge_const_handle aBisector ) const
  {
    CGAL_assertion(aBisector->is_bisector());
    CGAL_assertion(handle_assigned(aBisector->opposite()));
    CGAL_assertion(aBisector->opposite()->is_bisector());

    Halfedge_const_handle lBorderA = aBisector->defining_border();
    Halfedge_const_handle lBorderB = aBisector->opposite()->defining_border();

    return Construct_offset_point_2<Traits>(mTraits)()(aT,GetEdge(lBorderA),GetEdge(lBorderB));
  }

  Traits const&           mTraits ;
  boost::dynamic_bitset<> mVisitedBisectors;
  Halfedge_vector         mBorders ;
};

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Polygon_offset_builder_2.C>
#endif

#endif // CGAL_POLYGON_OFFSET_BUILDER_2_H //
// EOF //

