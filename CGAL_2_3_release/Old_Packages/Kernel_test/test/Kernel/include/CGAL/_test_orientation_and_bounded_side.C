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
// file          : _test_orientation_and_bounded_side.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_ORIENTATION_AND_BOUNDED_SIDE_C
#define CGAL__TEST_ORIENTATION_AND_BOUNDED_SIDE_C

#ifndef CGAL__TEST_ORIENTATION_AND_BOUNDED_SIDE_H
#include <CGAL/_test_orientation_and_bounded_side.h>
#endif // CGAL__TEST_ORIENTATION_AND_BOUNDED_SIDE_H


template <class R>
bool
_test_orientation_and_bounded_side(const R&)
{
  typedef typename R::RT    RT;
  typedef typename R::FT    FT;
  const RT RT0(0);
  const RT RT1(1);
  typedef CGAL::Point_2<R>  Point_2;
  typedef CGAL::Point_3<R>  Point_3;

  Point_2 e0( RT1, RT0);
  Point_2 e1( RT0, RT1);
  Point_2 e(-RT1, RT0);
  Point_2 org( RT0, RT0);
  assert( CGAL::orientation( org, e0, e1) == CGAL::POSITIVE );
  assert( CGAL::orientation( e0, e1, e) == CGAL::POSITIVE );
  assert( CGAL::side_of_oriented_circle( e0, e1, e, org) == \
            CGAL::ON_POSITIVE_SIDE );

  Point_3 s0( RT1, RT0, RT0);
  Point_3 s1( RT0, RT1, RT0);
  Point_3 s2( RT0, RT0, RT1);
  Point_3 s(-RT1, RT0, RT0);
  Point_3 org3( RT0, RT0, RT0);
  assert( CGAL::orientation( org3, s0, s1, s2) == CGAL::POSITIVE );
  assert( CGAL::orientation( s, s0, s1, s2) == CGAL::POSITIVE );
  assert( CGAL::side_of_oriented_sphere( s, s0, s1, s2, org3) == \
            CGAL::ON_POSITIVE_SIDE );

  return true;
}


#endif // CGAL__TEST_ORIENTATION_AND_BOUNDED_SIDE_C
