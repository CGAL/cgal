// Copyright (c) 2017
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL__TEST_FCT_WEIGHTED_POINT_2_H
#define CGAL__TEST_FCT_WEIGHTED_POINT_2_H

#include <CGAL/Line_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Weighted_point_2.h>
#include <CGAL/use.h>

#include <cassert>
#include <iostream>

#include "_approx_equal.h"

template <class R>
bool
_test_fct_weighted_point_2(const R& )
{
  std::cout << "Testing functions Weighted_point_2" ;

  typedef typename  R::RT    RT;

  CGAL::Point_2<R> p1( RT(18), RT(12), RT(3) );  // ( 6, 4)
  CGAL::Point_2<R> p2( RT(18), RT(15), RT(3) );  // ( 6, 5)
  CGAL::Point_2<R> p3( RT(18), RT( 9), RT(3) );  // ( 6, 3)
  CGAL::Point_2<R> p4( RT(28), RT(40), RT(4) );  // ( 7,10)
  CGAL::Point_2<R> p5( RT(12), RT(-4), RT(4) );  // ( 3,-1)
  CGAL::Point_2<R> p6( RT(28), RT(12), RT(4) );  // ( 7, 3)
  CGAL::Point_2<R> p7( RT(18), RT( 6), RT(3) );  // ( 6, 2)
  CGAL::Point_2<R> p8( RT(24), RT( 9), RT(3) );  // ( 8, 3)
  CGAL::Point_2<R> p9( RT( 6), RT(10), RT(1) );  // ( 6,10)
  CGAL::Point_2<R> p10( RT(-3), RT(8), RT(1) );  // ( -3,8)

  CGAL::Weighted_point_2<R> wp1( p1, 0);
  CGAL::Weighted_point_2<R> wp2( p2, 0);
  CGAL::Weighted_point_2<R> wp3( p3, 0);
  CGAL::Weighted_point_2<R> wp4( p4, 0);
  CGAL::Weighted_point_2<R> wp5( p5, 0);
  CGAL::Weighted_point_2<R> wp6( p6, 0);
  CGAL::Weighted_point_2<R> wp7( p7, 0);
  CGAL::Weighted_point_2<R> wp8( p8, 0);
  CGAL::Weighted_point_2<R> wp9( p9, 0);

  CGAL::Weighted_point_2<R> wp1_b( p1, 4);
  CGAL::Weighted_point_2<R> wp2_b( p2, -2);
  CGAL::Weighted_point_2<R> wp3_b( p3, 3);
  CGAL::Weighted_point_2<R> wp4_b( p4, 5);
  CGAL::Weighted_point_2<R> wp5_b( p5, 28);
  CGAL::Weighted_point_2<R> wp6_b( p6, -8);
  CGAL::Weighted_point_2<R> wp7_b( p7, 8);
  CGAL::Weighted_point_2<R> wp8_b( p8, 6);
  CGAL::Weighted_point_2<R> wp10_b( p10, -8);

  CGAL::Point_2<R> p_00(RT(0), RT(0) ); // 0, 0
  CGAL::Point_2<R> p_10(RT(4), RT(0), RT(1) ); // 4,0
  CGAL::Point_2<R> p_01(RT(0), RT(5), RT(1) ); // 0,5
  CGAL::Weighted_point_2<R> wp_00( p_00, RT(0) );
  CGAL::Weighted_point_2<R> wp_10( p_10, RT(16) );
  CGAL::Weighted_point_2<R> wp_01( p_01, RT(25) );

  assert( CGAL::compare_power_distance(p1, wp1, wp1) == CGAL::compare_distance(p1, p1, p1));
  assert( CGAL::compare_power_distance(p1, wp2, wp4) == CGAL::compare_distance(p1, p2, p4));

  assert( CGAL::compare_power_distance(p7, wp9, wp3) == CGAL::LARGER);
  assert( CGAL::compare_power_distance(p7, wp9, wp7) == CGAL::LARGER);
  assert( CGAL::compare_power_distance(p6, wp3, wp8) == CGAL::EQUAL);
  assert( CGAL::compare_power_distance(p2, wp5, wp5) == CGAL::EQUAL);
  assert( CGAL::compare_power_distance(p8, wp6, wp9) == CGAL::SMALLER);
  assert( CGAL::compare_power_distance(p6, wp5, wp4) == CGAL::SMALLER);

  assert( CGAL::compare_power_distance(p3, wp4_b, wp5_b) == CGAL::LARGER);
  assert( CGAL::compare_power_distance(p6, wp6_b, wp2_b) == CGAL::LARGER);
  assert( CGAL::compare_power_distance(p1, wp7_b, wp1_b) == CGAL::EQUAL);
  assert( CGAL::compare_power_distance(p3, wp3_b, wp5_b) == CGAL::EQUAL);
  assert( CGAL::compare_power_distance(p7, wp2_b, wp4_b) == CGAL::SMALLER);
  assert( CGAL::compare_power_distance(p2, wp1_b, wp2_b) == CGAL::SMALLER);

  std::cout << ".";

  assert( CGAL::power_product(wp1, wp1) == RT(0) );
  assert( CGAL::power_product(wp1, wp2) == CGAL::squared_distance(p1, p2) );
  assert( CGAL::power_product(wp3_b, wp5_b) == RT(-6) );

  std::cout << ".";

  CGAL::Line_2<R> l = CGAL::radical_axis(wp4_b, wp7_b);
  for(int i=-5; i<5; ++i)
  {
    CGAL::Point_2<R> p = R().construct_point_on_2_object()(l, i);
    assert(CGAL::compare_power_distance(p, wp4_b, wp7_b) == CGAL::EQUAL);
  }

  std::cout << ".";

  assert(power_side_of_oriented_power_circle(wp1, wp1_b) == CGAL::ON_POSITIVE_SIDE );
  assert(power_side_of_oriented_power_circle(wp2_b, wp2_b) == CGAL::ON_ORIENTED_BOUNDARY );
  assert(power_side_of_oriented_power_circle(wp6, wp6_b) == CGAL::ON_NEGATIVE_SIDE );

  assert(power_side_of_oriented_power_circle(wp7_b, wp4_b, wp8_b) == CGAL::ON_POSITIVE_SIDE );
  assert(power_side_of_oriented_power_circle(wp7_b, wp4_b, wp7_b) == CGAL::ON_ORIENTED_BOUNDARY );
  assert(power_side_of_oriented_power_circle(wp7_b, wp4_b, wp5_b) == CGAL::ON_NEGATIVE_SIDE );

  assert(power_side_of_oriented_power_circle(wp1_b, wp4_b, wp5_b, wp10_b) == CGAL::ON_POSITIVE_SIDE );
  assert(power_side_of_oriented_power_circle(wp1_b, wp4_b, wp5_b, wp5_b) == CGAL::ON_ORIENTED_BOUNDARY );
  assert(power_side_of_oriented_power_circle(wp1_b, wp4_b, wp5_b, wp6_b) == CGAL::ON_NEGATIVE_SIDE );

  std::cout << ".";

  assert( CGAL::power_side_of_bounded_power_circle(wp5_b, wp7_b) == CGAL::ON_UNBOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_circle(wp1_b, wp1_b) == CGAL::ON_BOUNDARY );
  assert( CGAL::power_side_of_bounded_power_circle(wp6_b, wp6) == CGAL::ON_BOUNDED_SIDE );

  assert( CGAL::power_side_of_bounded_power_circle(wp1_b, wp3_b, wp6_b) == CGAL::ON_UNBOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_circle(wp3_b, wp1_b, wp6_b) == CGAL::ON_UNBOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_circle(wp3_b, wp6_b, wp1_b) == CGAL::ON_BOUNDARY );
  assert( CGAL::power_side_of_bounded_power_circle(wp6_b, wp3_b, wp8_b) == CGAL::ON_BOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_circle(wp3_b, wp6_b, wp8_b) == CGAL::ON_BOUNDED_SIDE );

  assert( CGAL::power_side_of_bounded_power_circle(wp2, wp3, wp4, wp5)
            == CGAL::side_of_bounded_circle(p2, p3, p4, p5) );
  assert( CGAL::power_side_of_bounded_power_circle(wp2, wp3, wp4, wp4)
            == CGAL::side_of_bounded_circle(p2, p3, p4, p4) );

  assert( CGAL::power_side_of_bounded_power_circle(wp6_b, wp3_b, wp1_b, wp10_b) == CGAL::ON_UNBOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_circle(wp1_b, wp6_b, wp3_b, wp10_b) == CGAL::ON_UNBOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_circle(wp6_b, wp3_b, wp1_b, wp6_b) == CGAL::ON_BOUNDARY );
  assert( CGAL::power_side_of_bounded_power_circle(wp6_b, wp3_b, wp1_b, wp8_b) == CGAL::ON_BOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_circle(wp3_b, wp1_b, wp6_b, wp8_b) == CGAL::ON_BOUNDED_SIDE );

  std::cout << ".";

  assert( CGAL::squared_radius_smallest_orthogonal_circle(wp1) == RT(0) );
  assert( CGAL::squared_radius_smallest_orthogonal_circle(wp1_b) == -wp1_b.weight() );

  std::cout << CGAL::squared_radius_smallest_orthogonal_circle(wp1_b, wp3_b) << std::endl;

  assert( CGAL::squared_radius_smallest_orthogonal_circle(wp3, wp8) == RT(1) );
  assert( CGAL::squared_radius_smallest_orthogonal_circle(wp1_b, wp3_b) == RT(-3));

  std::cout << CGAL::weighted_circumcenter(wp_00, wp_10, wp_01) << std::endl;

  using CGAL::testsuite::approx_equal;
  assert( approx_equal(CGAL::squared_radius_smallest_orthogonal_circle(wp1, wp3, wp5),
                       CGAL::squared_radius(p1, p3, p5)) );
  assert( CGAL::squared_radius_smallest_orthogonal_circle(wp_00, wp_10, wp_01) == RT(0));

  std::cout << "done" << std::endl;
  return true;
}

#endif // CGAL__TEST_FCT_WEIGHTED_POINT_2_H
