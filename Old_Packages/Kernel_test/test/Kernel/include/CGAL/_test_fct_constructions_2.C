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
// file          : _test_fct_constructions_2.C
// revision      : 3.8
// revision_date : 08 Oct 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_FCT_CONSTRUCTIONS_2_C
#define CGAL__TEST_FCT_CONSTRUCTIONS_2_C

#ifndef CGAL__TEST_FCT_CONSTRUCTIONS_2_H
#include <CGAL/_test_fct_constructions_2.h>
#endif // CGAL__TEST_FCT_CONSTRUCTIONS_2_H

template <class R>
bool
_test_fct_constructions_2(const R&)
{
  typedef typename R::RT     RT;
  typedef CGAL::Point_2<R>    Point;
  typedef CGAL::Vector_2<R>   Vector;

  RT RT0(0);
  RT RT1(1);
  RT RT2(2);
  RT RT3(3);
  RT RT4(4);
  RT RT8(8);

  Point p( RT4, RT8, RT2);   // ( 2, 4)
  Point pne = p + Vector( RT1, RT1 );
  Point pnw = p + Vector(-RT1, RT1 );
  Point pse = p + Vector( RT1,-RT1 );
  Point psw = p + Vector(-RT1,-RT1 );

  // midpoint
  assert( CGAL::midpoint( pne, psw) == p);
  assert( CGAL::midpoint( pnw, pse) == p);

  // circumcenter
  assert( CGAL::circumcenter( pne, pse, pnw) == p);
  assert( CGAL::circumcenter( psw, pse, pnw) == p);

  // general position intersection point

  return true;
}


#endif // CGAL__TEST_FCT_CONSTRUCTIONS_2_C
