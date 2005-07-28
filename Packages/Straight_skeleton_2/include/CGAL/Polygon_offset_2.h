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
// file          : include/CGAL/Polygon_offset_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================

#ifndef CGAL_POLYGON_OFFSET_2_H
#define CGAL_POLYGON_OFFSET_2_H 1

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/optional.hpp>

#ifndef CGAL_NUMBER_UTILS_CLASSES_H
#include <CGAL/Number_utils_classes.h>
#endif

#ifndef CGAL_POLYGON_OFFSET_TRAITS_2_H
#include <CGAL/Polygon_offset_traits_2.h>
#endif

CGAL_BEGIN_NAMESPACE

template<class OffsetGraph_, class Traits_>
class Polygon_offset_2 
{
public :

  typedef OffsetGraph_ Graph ;
  typedef Traits_      Traits ;
 
  typedef typename Traits::FT FT ;
  
  typedef typename Traits::Polygon_2 Polygon_2 ;
  
  typedef typename Polygon_2::Segment_2 Segment_2 ;
  typedef typename Polygon_2::Point_2   Point_2 ;

  typedef boost::shared_ptr<Polygon_2> Polygon_2_Ptr ;
    
  Polygon_offset_2( Graph const& aGraph, Traits const& aTraits = Traits() )  ;

  template<class OutputIterator>
  OutputIterator Create ( FT aTime, OutputIterator aOut ) ;
  
private:

  typedef typename Graph::Halfedge_const_handle Halfedge_const_handle  ;
  
  typedef std::vector<Halfedge_const_handle> Halfedge_vector ;

  typedef boost::optional<Point_2> Optional_point_2 ;
    
  Min<FT> Min ;
  Max<FT> Max ;
  
  bool handled_assigned( Halfedge_const_handle aH ) const
  { 
    const Halfedge_const_handle cNull ;
    return aH != cNull ;
  }
  
  Halfedge_const_handle LocateHook( FT aTime, Halfedge_const_handle aBisector ) ;
  
  template<class OutputIterator>
  OutputIterator TraceOffsetPolygon( FT aTime, Halfedge_const_handle aHook, OutputIterator aOut ) ;
  
  Halfedge_const_handle LocateSeed( FT aTime ) ;
  
  bool IsVisited( Halfedge_const_handle aBisector ) { return mVisitedBisectors[aBisector->id()]; }
  
  void Visit( Halfedge_const_handle aBisector ) { mVisitedBisectors[aBisector->id()] = true ; }
  
  bool IsTimeWithinBisector ( FT aTime, Halfedge_const_handle aBisector )
  {
    CGAL_assertion(aBisector->is_bisector());
    CGAL_assertion(handle_assigned(aBisector->opposite()));
    CGAL_assertion(aBisector->opposite()->is_bisector());
    
    FT lTimeA = aBisector->vertex()->time();
    FT lTimeB = aBisector->opposite()->vertex()->time();
  
    FT lLo = Min(lTimeA,lTimeB);
    FT lHi = Max(lTimeA,lTimeB);
  
    return lLo < aTime && lHi >= aTime ;
  }
  
  Traits const&           mTraits ;
  boost::dynamic_bitset<> mVisitedBisectors;
  Halfedge_vector         mBorders ;
};

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Polygon_offset_2.C>
#endif

#endif // CGAL_POLYGON_OFFSET_2_H //
// EOF //
 
