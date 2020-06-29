// Copyright (c) 2001
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL__TEST_ANGLE_H
#define CGAL__TEST_ANGLE_H

template <class R>
bool
_test_angle(const R&)
{
  typedef typename R::RT    RT;
  const RT RT0(0);
  const RT RT1(1);
  typedef CGAL::Point_2<R>  Point_2;
  typedef CGAL::Point_3<R>  Point_3;

  typedef CGAL::Vector_2<R> Vector_2;
  typedef CGAL::Vector_3<R> Vector_3;

  Point_2 p(RT(2),RT(1));
  Point_2 q(RT(5),RT(4));
  Point_2 r(RT(5),RT(10));
  Point_2 s(RT0, RT0);

  Vector_2 qp = p - q;
  Vector_2 qr = r - q;

  assert( CGAL::angle( p, q, r ) == CGAL::OBTUSE );
  assert( CGAL::angle( qp , qr ) == CGAL::OBTUSE );
  assert( CGAL::angle( q, p, q, r ) == CGAL::OBTUSE );
  assert( CGAL::angle( r, p, q ) == CGAL::ACUTE );
  assert( CGAL::angle( p, s, q , r) == CGAL::OBTUSE );
  assert( CGAL::angle( p, r, s , q) == CGAL::ACUTE );

  Point_2 e0( RT1, RT0);
  Point_2 e1( RT0, RT1);
  Point_2 org( RT0, RT0);
  assert( CGAL::angle( e0, org, e1) == CGAL::RIGHT );

  Point_3 sx( RT1, RT0, RT0);
  Point_3 msx( - RT1, RT0, RT0);
  Point_3 sy( RT0, RT1, RT0);
  Point_3 sz( RT0, RT0, RT1);
  Point_3 org3( RT0, RT0, RT0);
  assert( CGAL::angle( sx, org3, sy ) == CGAL::RIGHT );
  assert( CGAL::angle( sx, sy, sz) == CGAL::ACUTE );
  assert( CGAL::angle( org3, sy, sx, sz) == CGAL::RIGHT );
  assert( CGAL::angle( org3, sz, sz, sy) == CGAL::OBTUSE );
  assert( CGAL::angle( org3, sx, sy, sx) == CGAL::ACUTE );
  assert( CGAL::angle( org3 - sy, sx - sz) == CGAL::RIGHT );
  assert( CGAL::angle( org3 - sz, sz - sy) == CGAL::OBTUSE );
  assert( CGAL::angle( org3 - sx, sy - sx) == CGAL::ACUTE );

  Vector_3 vz(RT0, RT0, RT1);
  Vector_3 vcoplanar(RT1, RT1, RT0);
  Vector_3 vmz(RT0, RT0, -RT1);
  assert( CGAL::angle( org3, sx, sy, vz) ==  CGAL::ACUTE );
  assert( CGAL::angle( org3, sx, sy, vcoplanar) ==  CGAL::RIGHT );
  assert( CGAL::angle( org3, sx, sy, vmz) ==  CGAL::OBTUSE );
  assert( CGAL::approximate_angle( sx, org3, sy ) > RT(89) );
  assert( CGAL::approximate_angle( sx-org3, sy-org3 ) > RT(89) );
  assert( CGAL::approximate_angle( sx, org3, sy ) < RT(91) );
  assert( CGAL::approximate_angle( sx, org3, sx ) == RT(0) );
  assert( CGAL::approximate_angle( sx, org3, msx ) > RT(179) );
  assert( CGAL::approximate_angle( sx, org3, msx ) < RT(181) );
  return true;
}

#endif // CGAL__TEST_ANGLE_H
