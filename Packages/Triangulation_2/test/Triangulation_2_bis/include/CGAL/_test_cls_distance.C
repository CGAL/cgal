// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL
// release
// of the Computational Geometry Algorithms Library (CGAL). It is
// not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : 
// file          : 
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat
// (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <cassert>

template<class Traits, class Point>
void
_test_cls_distance(Point p[34],const Traits &)
{

  typedef typename Traits::Distance Distance;
  //typedef typename Traits::Point  Point;
  //test for constructors
  Distance d0;
  Distance d1(p[0]); 
  Distance d2(p[0],p[1]); 
  Distance d3(p[0],p[1],p[2]);
  // try to apply compare to a valid objects
   assert(d3.compare()==CGAL::SMALLER);
   // test accessors
   assert(d1.get_point(0)==p[0]);
   assert(d3.get_point(0)==p[0]);
   assert(d3.get_point(1)==p[1]);
   assert(d3.get_point(2)==p[2]);
   // test modifiers
   // Not possible to apply set_point to an empty Distance.
   d0.set_point(0,p[0]);
   d0.set_point(1,p[1]);
   d0.set_point(2,p[2]);
   assert(d0.compare()==CGAL::SMALLER);
   d1.set_point(1,p[1]);
   d1.set_point(2,p[2]);
   assert(d1.compare()==CGAL::SMALLER);
   // test a degenerate case p1=p2
   d3.set_point(2,p[1]);
   assert(d3.compare()==CGAL::EQUAL);
   // test a degenerate case p0=p1=p2
   d3.set_point(0,p[1]);
   assert(d3.compare()==CGAL::EQUAL);
 }
