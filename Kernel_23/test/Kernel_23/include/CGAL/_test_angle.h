// Copyright (c) 2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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
  typedef typename R::FT    FT;
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

  Point_3 s0( RT1, RT0, RT0);
  Point_3 s1( RT0, RT1, RT0);
  Point_3 s2( RT0, RT0, RT1);
  Point_3 org3( RT0, RT0, RT0);
  assert( CGAL::angle( s0, org3, s1 ) == CGAL::RIGHT );
  assert( CGAL::angle( s0, s1, s2) == CGAL::ACUTE );

  return true;
}

#endif // CGAL__TEST_ANGLE_H
