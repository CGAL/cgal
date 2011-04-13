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
// file          : _test_fct_constructions_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_FCT_CONSTRUCTIONS_3_C
#define CGAL__TEST_FCT_CONSTRUCTIONS_3_C

#include <CGAL/_test_fct_constructions_3.h>

template <class R>
bool
_test_fct_constructions_3(const R&)
{
  typedef typename R::RT     RT;
  typedef typename R::Point_3   Point;
  typedef typename R::Vector_3  Vector;

  RT RT0(0);
  RT RT1(1);
  RT RT2(2);
  RT RT3(3);
  RT RT4(4);
  RT RT8(8);

  Point p( RT4, RT8, -RT2, RT2);   // ( 2, 4, -1)
  Point p111 = p + Vector( RT1, RT1, RT1 );
  Point p011 = p + Vector(-RT1, RT1, RT1 );
  Point p101 = p + Vector( RT1,-RT1, RT1 );
  Point p001 = p + Vector(-RT1,-RT1, RT1 );
  Point p000 = p + Vector(-RT1,-RT1,-RT1 );
  Point p100 = p + Vector( RT1,-RT1,-RT1 );
  Point p110 = p + Vector( RT1, RT1,-RT1 );
  Point p010 = p + Vector(-RT1, RT1,-RT1 );

  Point p2   = p + Vector(-RT1, RT0, RT0 );
  Point p3   = p + Vector( RT1, RT0, RT0 );
  Point p4   = p + Vector( RT0, RT1, RT0 );

  // midpoint
  assert( CGAL::midpoint( p111, p000) == p);
  assert( CGAL::midpoint( p110, p001) == p);
  assert( CGAL::midpoint( p010, p101) == p);
  assert( CGAL::midpoint( p100, p011) == p);

  // circumcenter
  assert( CGAL::circumcenter( p111, p001, p010, p000) == p);
  assert( CGAL::circumcenter( p101, p001, p010, p100) == p);
  assert( CGAL::circumcenter( p001, p000, p110, p100) == p);

  assert( CGAL::circumcenter( p2, p3, p4) == p);

  // centroid
  Point p_11 = p + Vector(RT0, RT1, RT1);
  assert( CGAL::centroid( p111, p010, p101, p000) == p);
  assert( CGAL::centroid( p111, p_11, p011 ) == p_11);

  // projection onto a plane

  return true;
}

#endif // CGAL__TEST_FCT_CONSTRUCTIONS_3_C
