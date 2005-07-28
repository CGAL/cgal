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
// file          : include/CGAL/POlygon_offset_traits_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_POLYGON_OFFSET_TRAITS_2_H
#define CGAL_POLYGON_OFFSET_TRAITS_2_H 1

#include <CGAL/Polygon_2.h>

#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>

CGAL_BEGIN_NAMESPACE

template<class R>
class Polygon_offset_traits_2 
{
public:
  
  typedef R Rep ;

  typedef Polygon_2<R> Polygon_2 ;
  
  typedef typename Polygon_2::Point_2 Point_2 ;
  typedef typename Polygon_2::FT      FT ;

  typedef boost::shared_ptr<Polygon_2> Polygon_2_ptr ;  
  
  typedef boost::optional<Point_2> Optional_point_2 ;
  
public:

  
  Polygon_2_ptr Construct_polygon_2_ptr() const { return Polygon_2_ptr(new Polygon_2()); }
  
  Optional_point_2 Construct_offset_point_2( FT             aTime
                                           , Point_2 const& aBorderS
                                           , Point_2 const& aBorderT
                                           , Point_2 const& aBisectorS
                                           , Point_2 const& aBisectorT
                                           ) const ; 
   
private :

  typedef typename Rep::Line_2  Line_2 ;
  
  typename Rep::Construct_line_2 construct_line ;

  Line_2 Construct_offset_line( Line_2 const& aLine, FT aSquaredOffsetDist ) const ;
  
  Optional_point_2 Intersect( Line_2 const& aA, Line_2 const& aB ) const ;
                             
} ;

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Polygon_offset_traits_2.C>
#endif


#endif // CGAL_POLYGON_OFFSET_TRAITS_2_H //
// EOF //
