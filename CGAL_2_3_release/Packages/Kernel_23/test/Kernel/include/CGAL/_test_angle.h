// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : _test_angle.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 
#ifndef CGAL__TEST_ANGLE_H
#define CGAL__TEST_ANGLE_H

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/predicates_on_points_3.h>

template <class R>
bool
_test_angle(const R&)
{
  typedef typename R::RT    RT;
  typedef typename R::FT    FT;
  const RT RT0(0);
  const RT RT1(1);
  typedef CGAL::Point_2<R>  Point_2;
  typedef CGAL::Point_3<R>  Point_3;

  Point_2 p(RT(2),RT(1));
  Point_2 q(RT(5),RT(4));
  Point_2 r(RT(5),RT(10));

  assert( CGAL::angle( p, q, r ) == CGAL::OBTUSE );
  assert( CGAL::angle( r, p, q ) == CGAL::ACUTE );

  Point_2 e0( RT1, RT0);
  Point_2 e1( RT0, RT1);
  Point_2 org( RT0, RT0);
  assert( CGAL::angle( e0, org, e1) == CGAL::RIGHT );

  Point_3 s0( RT1, RT0, RT0);
  Point_3 s1( RT0, RT1, RT0);
  Point_3 s2( RT0, RT0, RT1);
  Point_3 org3( RT0, RT0, RT0);
  assert( CGAL::angle( s0, org3, s1 ) == CGAL::RIGHT );
  assert( CGAL::angle( s0, s1, s2) == CGAL::ACUTE );

  return true;
}

#endif // CGAL__TEST_ANGLE_H
