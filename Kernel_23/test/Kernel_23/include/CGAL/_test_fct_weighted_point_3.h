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

#ifndef CGAL__TEST_FCT_WEIGHTED_POINT_3_H
#define CGAL__TEST_FCT_WEIGHTED_POINT_3_H

#include <CGAL/Segment_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Weighted_point_3.h>

#include <cassert>
#include <iostream>

template <class R>
bool
_test_fct_weighted_point_3(const R& )
{
  std::cout << "Testing functions Weighted_point_3" ;

  typedef typename  R::RT    RT;
  const bool nonexact = std::is_floating_point<RT>::value;

  CGAL::Point_3<R> p1(RT(18), RT(15), RT(-21), RT(3) ); //  6,  5, -7
  CGAL::Point_3<R> p2(RT(18), RT(15), RT( 12), RT(3) ); //  6,  5,  4
  CGAL::Point_3<R> p3(RT(18), RT(12), RT(-21), RT(3) ); //  6,  4, -7
  CGAL::Point_3<R> p4(RT(28), RT(40), RT( 20), RT(4) ); //  7, 10,  5
  CGAL::Point_3<R> p5(RT(12), RT(-4), RT(-20), RT(4) ); //  3, -1, -5
  CGAL::Point_3<R> p6(RT(18), RT(18), RT(-21), RT(3) ); //  6,  6, -7
  CGAL::Point_3<R> p7(RT(16), RT(15),  RT(-7), RT(3) ); // barycenter of p4 p5 p6

  CGAL::Weighted_point_3<R> wp1( p1, RT(0));
  CGAL::Weighted_point_3<R> wp2( p2, RT(0));
  CGAL::Weighted_point_3<R> wp3( p3, RT(0));
  CGAL::Weighted_point_3<R> wp4( p4, RT(0));
  CGAL::Weighted_point_3<R> wp5( p5, RT(0));
  CGAL::Weighted_point_3<R> wp6( p6, RT(0));
  CGAL::Weighted_point_3<R> wp7( p7, RT(0));

  CGAL::Weighted_point_3<R> wp1_b( p1, RT(-20));
  CGAL::Weighted_point_3<R> wp2_b( p2, RT(-10));
  CGAL::Weighted_point_3<R> wp3_b( p3, RT(5));
  CGAL::Weighted_point_3<R> wp4_b( p4, RT(174));
  CGAL::Weighted_point_3<R> wp5_b( p5, RT(-3));
  CGAL::Weighted_point_3<R> wp6_b( p6, RT(0));
  CGAL::Weighted_point_3<R> wp6b_b( p6, RT(4));

  CGAL::Point_3<R> p000(RT(0), RT(0), RT(0) ); // 0, 0, 0
  CGAL::Point_3<R> p100(RT(4), RT(0), RT(0), RT(1) ); // 4,0,0
  CGAL::Point_3<R> p010(RT(0), RT(5), RT(0), RT(1) ); // 0,5,0
  CGAL::Point_3<R> p001(RT(0), RT(0), RT(-6), RT(1) ); // 0,0,6
  CGAL::Point_3<R> p011(RT(0), RT(3), RT(2), RT(1) ); // 0,3,2
  CGAL::Weighted_point_3<R> wp000( p000, RT(0) );
  CGAL::Weighted_point_3<R> wp100( p100, RT(16) );
  CGAL::Weighted_point_3<R> wp010( p010, RT(25) );
  CGAL::Weighted_point_3<R> wp001( p001, RT(36) );
  CGAL::Weighted_point_3<R> wp011( p011, RT(14) );
  CGAL::Weighted_point_3<R> wp000m( p000, RT(100) );

  assert( CGAL::compare_power_distance(p1, wp1, wp1) == CGAL::compare_distance(p1, p1, p1));
  assert( CGAL::compare_power_distance(p1, wp2, wp4) == CGAL::compare_distance(p1, p2, p4));

  assert( CGAL::compare_power_distance(p3, wp6, wp5) == CGAL::SMALLER);
  assert( CGAL::compare_power_distance(p5, wp5, wp3) == CGAL::SMALLER);
  assert( CGAL::compare_power_distance(p1, wp3, wp6) == CGAL::EQUAL);
  assert( CGAL::compare_power_distance(p2, wp2, wp2) == CGAL::EQUAL);
  assert( CGAL::compare_power_distance(p2, wp5, wp1) == CGAL::LARGER);
  assert( CGAL::compare_power_distance(p4, wp2, wp4) == CGAL::LARGER);

  assert( CGAL::compare_power_distance(p2, wp2_b, wp4  ) == CGAL::SMALLER);
  assert( CGAL::compare_power_distance(p4, wp3_b, wp1_b) == CGAL::SMALLER);
  assert( CGAL::compare_power_distance(p2, wp3_b, wp3_b) == CGAL::EQUAL);
  assert( CGAL::compare_power_distance(p1, wp4_b, wp3_b) == CGAL::EQUAL);
  assert( CGAL::compare_power_distance(p3, wp5_b, wp1_b) == CGAL::LARGER);
  assert( CGAL::compare_power_distance(p5, wp2_b, wp6_b) == CGAL::LARGER);

  std::cout << ".";

  assert( CGAL::power_product(wp1, wp1) == RT(0) );
  assert( CGAL::power_product(wp1, wp2) == CGAL::squared_distance(p1, p2) );
  assert( CGAL::power_product(wp3_b, wp5_b) == RT(36) );

  std::cout << ".";

  assert( CGAL::power_side_of_oriented_power_sphere(wp6b_b, wp6_b) == CGAL::ON_NEGATIVE_SIDE );
  assert( CGAL::power_side_of_oriented_power_sphere(wp2_b, wp2_b) == CGAL::ON_ORIENTED_BOUNDARY );
  assert( CGAL::power_side_of_oriented_power_sphere(wp6_b, wp6b_b) == CGAL::ON_POSITIVE_SIDE );

  assert( CGAL::power_side_of_oriented_power_sphere(wp3, wp1, wp6) == CGAL::ON_NEGATIVE_SIDE
          && ! CGAL::Segment_3<R>(p3, p1).has_on(p6) );
  assert( CGAL::power_side_of_oriented_power_sphere(wp5, wp6, wp5) == CGAL::ON_ORIENTED_BOUNDARY
          && CGAL::Segment_3<R>(p5, p6).has_on(p5) );
  assert( CGAL::power_side_of_oriented_power_sphere(wp3, wp6, wp1) == CGAL::ON_POSITIVE_SIDE
          && CGAL::Segment_3<R>(p3, p6).has_on(p1) );

  assert( CGAL::power_side_of_oriented_power_sphere(wp3_b, wp6b_b, wp1_b) == CGAL::ON_NEGATIVE_SIDE );
  assert( CGAL::power_side_of_oriented_power_sphere(wp1_b, wp3_b, wp3_b) == CGAL::ON_ORIENTED_BOUNDARY );
  assert( CGAL::power_side_of_oriented_power_sphere(wp3_b, wp1_b, wp6b_b) == CGAL::ON_POSITIVE_SIDE );

  assert( CGAL::power_side_of_oriented_power_sphere(wp000, wp100, wp011, wp010) == CGAL::ON_NEGATIVE_SIDE );
  assert( CGAL::power_side_of_oriented_power_sphere(wp000, wp100, wp010, wp100) == CGAL::ON_ORIENTED_BOUNDARY );
  assert( CGAL::power_side_of_oriented_power_sphere(wp000, wp100, wp010, wp011) == CGAL::ON_POSITIVE_SIDE );

  assert( CGAL::power_side_of_oriented_power_sphere(wp2, wp3, wp4, wp5, wp6)
            == CGAL::side_of_oriented_sphere(p2, p3, p4, p5, p6));
  assert( CGAL::power_side_of_oriented_power_sphere(wp2, wp3, wp4, wp5, wp5)
            == CGAL::side_of_oriented_sphere(p2, p3, p4, p5, p5));

  assert( CGAL::power_side_of_oriented_power_sphere(wp000, wp100, wp010, wp001, wp000m) == CGAL::ON_NEGATIVE_SIDE );
  assert( CGAL::power_side_of_oriented_power_sphere(wp000, wp100, wp010, wp001, wp010) == CGAL::ON_ORIENTED_BOUNDARY );
  assert( CGAL::power_side_of_oriented_power_sphere(wp000, wp100, wp010, wp001, wp1_b) == CGAL::ON_POSITIVE_SIDE );

  std::cout << ".";

  assert( CGAL::power_side_of_bounded_power_sphere(wp6b_b, wp6_b) == CGAL::ON_UNBOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_sphere(wp2_b, wp2_b) == CGAL::ON_BOUNDARY );
  assert( CGAL::power_side_of_bounded_power_sphere(wp6_b, wp6b_b) == CGAL::ON_BOUNDED_SIDE );

  assert( CGAL::power_side_of_bounded_power_sphere(wp3_b, wp6b_b, wp1_b) == CGAL::ON_UNBOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_sphere(wp1_b, wp3_b, wp3_b) == CGAL::ON_BOUNDARY );
  assert( CGAL::power_side_of_bounded_power_sphere(wp3_b, wp1_b, wp6b_b) == CGAL::ON_BOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_sphere(wp1_b, wp3_b, wp6b_b) == CGAL::ON_BOUNDED_SIDE );

  assert( CGAL::coplanar(p4, p5, p6, p7) || nonexact );
  assert( CGAL::power_side_of_bounded_power_sphere(wp4, wp5, wp6, wp7)
          == CGAL::side_of_bounded_sphere(p4, p5, p6, p7) );
  assert( CGAL::power_side_of_bounded_power_sphere(wp4, wp5, wp6, wp5)
          == CGAL::side_of_bounded_sphere(p4, p5, p6, p5) );
  assert( CGAL::power_side_of_bounded_power_sphere(wp4, wp5, wp7, wp6)
          == CGAL::side_of_bounded_sphere(p4, p5, p7, p6) );

  assert( CGAL::power_side_of_bounded_power_sphere(wp000, wp100, wp011, wp010) == CGAL::ON_UNBOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_sphere(wp000, wp100, wp010, wp100) == CGAL::ON_BOUNDARY );
  assert( CGAL::power_side_of_bounded_power_sphere(wp000, wp100, wp010, wp011) == CGAL::ON_BOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_sphere(wp100, wp000, wp010, wp011) == CGAL::ON_BOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_sphere(wp010, wp000, wp100, wp011) == CGAL::ON_BOUNDED_SIDE );

  assert( CGAL::power_side_of_bounded_power_sphere(wp2, wp3, wp4, wp5, wp6)
            == CGAL::side_of_bounded_sphere(p2, p3, p4, p5, p6));
  assert( CGAL::power_side_of_bounded_power_sphere(wp2, wp3, wp4, wp5, wp5)
            == CGAL::side_of_bounded_sphere(p2, p3, p4, p5, p5));

  assert( CGAL::power_side_of_bounded_power_sphere(wp000, wp100, wp010, wp001, wp000m) == CGAL::ON_BOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_sphere(wp000, wp100, wp010, wp001, wp010) == CGAL::ON_BOUNDARY );
  assert( CGAL::power_side_of_bounded_power_sphere(wp000, wp100, wp010, wp001, wp1_b) == CGAL::ON_UNBOUNDED_SIDE );
  assert( CGAL::power_side_of_bounded_power_sphere(wp100, wp000, wp001, wp010, wp1_b) == CGAL::ON_UNBOUNDED_SIDE );

  std::cout << ".";

  assert( CGAL::squared_radius_smallest_orthogonal_sphere(wp1) == RT(0) );
  assert( CGAL::squared_radius_smallest_orthogonal_sphere(wp1_b) == -wp1_b.weight() );

  assert( CGAL::squared_radius_smallest_orthogonal_sphere(wp3, wp6) == RT(1) );
  assert( CGAL::squared_radius_smallest_orthogonal_sphere(wp1_b, wp3_b) == RT(164));

  assert( CGAL::squared_radius_smallest_orthogonal_sphere(wp1, wp3, wp5)
            == CGAL::squared_radius(p1, p3, p5) || nonexact);
  assert( CGAL::squared_radius_smallest_orthogonal_sphere(wp000, wp100, wp010) == RT(0));

  assert( CGAL::squared_radius_smallest_orthogonal_sphere(wp1, wp3, wp4, wp5)
            == CGAL::squared_radius(p1, p3, p4, p5) || nonexact);
  assert( CGAL::squared_radius_smallest_orthogonal_sphere(wp000, wp100, wp010, wp001) == RT(0));

  std::cout << ".";

  assert( CGAL::power_distance_to_power_sphere(wp1, wp3_b, wp5, wp4_b, wp1) == RT(0) );

  std::cout << ".";

  assert( CGAL::compare_weighted_squared_radius(wp4_b, RT(10)) == CGAL::SMALLER );
  assert( CGAL::compare_weighted_squared_radius(wp1, RT(0)) == CGAL::EQUAL );
  assert( CGAL::compare_weighted_squared_radius(wp1_b, RT(0)) == CGAL::LARGER );

  assert( CGAL::compare_weighted_squared_radius(wp1_b, wp2, RT(0)) == CGAL::LARGER );
  assert( CGAL::compare_weighted_squared_radius(wp3, wp6, RT(1)) == CGAL::EQUAL );
  assert( CGAL::compare_weighted_squared_radius(wp2, wp4_b, RT(255)) == CGAL::SMALLER );

  assert( CGAL::compare_weighted_squared_radius(wp1, wp2, wp3, RT(50)) == CGAL::SMALLER);
  assert( CGAL::compare_weighted_squared_radius(wp000, wp100, wp010, RT(0)) == CGAL::EQUAL);
  assert( CGAL::compare_weighted_squared_radius(wp4_b, wp2, wp3, RT(50)) == CGAL::LARGER );

  assert( CGAL::compare_weighted_squared_radius(wp000, wp100, wp010, wp2, RT(100)) == CGAL::SMALLER );
  assert( CGAL::compare_weighted_squared_radius(wp000, wp100, wp010, wp001, RT(0)) == CGAL::EQUAL );
  assert( CGAL::compare_weighted_squared_radius(wp000, wp100, wp010, wp3, RT(50)) == CGAL::LARGER );

  std::cout << "done" << std::endl;
  return true;
}

#endif // CGAL__TEST_FCT_WEIGHTED_POINT_3_H
