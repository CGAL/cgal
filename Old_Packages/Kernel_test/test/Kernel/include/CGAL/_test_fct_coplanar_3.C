// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        :
// file          : _test_fct_coplanar_3.C
// revision      : 3.8
// revision_date : 08 Oct 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_FCT_COPLANAR_3_C
#define CGAL__TEST_FCT_COPLANAR_3_C

#ifndef CGAL__TEST_FCT_COPLANAR_3_H
#include <CGAL/_test_fct_coplanar_3.h>
#endif // CGAL__TEST_FCT_COPLANAR_3_H

template <class R>
bool
_test_fct_coplanar_3(const R& )
{
  typedef typename R::RT     RT;
  typedef CGAL::Point_3<R>   Point;
  typedef CGAL::Vector_3<R>  Vector;
    RT RT0(0);
    RT RT1(1);
    RT RT2(2);
    RT RT3(3);
    RT RT4(4);
    RT RT6(6);
    RT RT8(8);
  
  Point p = Point( RT1, RT0, RT1, RT2);
  Point q = Point( RT4, RT1, RT2, RT8);
  Point r = Point( RT3, RT1, RT3, RT6);
  Point s = p + (q - r);
  assert( CGAL::coplanar( p,q,r,s));
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::NEGATIVE);
  assert( CGAL::coplanar_orientation( p,q,s,r) == CGAL::NEGATIVE );
  assert( CGAL::coplanar_orientation( p,q,r,r) == CGAL::POSITIVE );
  s = p + RT2*( q - p);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::COLLINEAR );
  s = p - (q - r);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::POSITIVE );
  p = Point( RT0, RT1, RT1, RT2);
  q = Point( RT1, RT4, RT2, RT8);
  r = Point( RT1, RT3, RT3, RT6);
  s = p + (q - r);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::NEGATIVE );
  assert( CGAL::coplanar_orientation( p,q,s,r) == CGAL::NEGATIVE );
  assert( CGAL::coplanar_orientation( p,q,r,r) == CGAL::POSITIVE );
  s = p + RT2*( q - p);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::COLLINEAR );
  s = p - (q - r);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::POSITIVE );
  p = Point( RT0, RT1, RT1, RT2);
  q = Point( RT1, RT2, RT4, RT8);
  r = Point( RT1, RT3, RT3, RT6);
  s = p + (q - r);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::NEGATIVE );
  assert( CGAL::coplanar_orientation( p,q,s,r) == CGAL::NEGATIVE );
  assert( CGAL::coplanar_orientation( p,q,r,r) == CGAL::POSITIVE );
  s = p + RT2*( q - p);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::COLLINEAR );
  s = p - (q - r);
  assert( CGAL::coplanar_orientation( p,q,r,s) == CGAL::POSITIVE );
  return true;
}


#endif // CGAL__TEST_FCT_COPLANAR_3_C
