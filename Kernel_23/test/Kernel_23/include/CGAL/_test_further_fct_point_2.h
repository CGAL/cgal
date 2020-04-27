// Copyright (c) 1999
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
// Author(s)     : Stefan Schirra


#ifndef CGAL__TEST_FURTHER_FCT_POINT_2_H
#define CGAL__TEST_FURTHER_FCT_POINT_2_H

#include "_approx_equal.h"
#include <boost/type_traits/is_same.hpp>

template <class R>
bool
_test_further_fct_point_2(const R& )
{
 std::cout << "Testing further functions Point_2" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 CGAL::Point_2<R> p1( RT(18), RT(12), RT(3) );  // ( 6, 4)
 CGAL::Point_2<R> p2( RT(18), RT(15), RT(3) );  // ( 6, 5)
 CGAL::Point_2<R> p3( RT(18), RT( 9), RT(3) );  // ( 6, 3)
 CGAL::Point_2<R> p4( RT(28), RT(40), RT(4) );  // ( 7,10)
 CGAL::Point_2<R> p5( RT(12), RT(-4), RT(4) );  // ( 3,-1)
 CGAL::Point_2<R> p6( RT(28), RT(12), RT(4) );  // ( 7, 3)
 CGAL::Point_2<R> p7( RT(18), RT( 6), RT(3) );  // ( 6, 2)
 CGAL::Point_2<R> p8( RT(24), RT( 9), RT(3) );  // ( 8, 3)
 CGAL::Point_2<R> p9( RT( 6), RT(10), RT(1) );  // ( 6,10)

#if 0
 // This function is not documented, not used in CGAL, and has no functor.
 assert( CGAL::compare_deltax_deltay( p1, p4, p5, p6) == CGAL::SMALLER);
 assert( CGAL::compare_deltax_deltay( p5, p6, p2, p1) == CGAL::LARGER);
 assert( CGAL::compare_deltax_deltay( p6, p5, p2, p1) == CGAL::LARGER);
 assert( CGAL::compare_deltax_deltay( p1, p2, p5, p6) == CGAL::SMALLER);
 assert( CGAL::compare_deltax_deltay( p4, p3, p3, p7) == CGAL::EQUAL);
 assert( CGAL::compare_deltax_deltay( p2, p3, p3, p6) == CGAL::EQUAL);
#endif

 std::cout << '.';

 assert( CGAL::x_equal(p1,p1) );
 assert( CGAL::x_equal(p2,p3) );
 assert( !CGAL::x_equal(p2,p4) );

 assert( CGAL::y_equal(p3,p6) );
 assert( !CGAL::y_equal(p1,p3) );

 assert( CGAL::compare_x( p1, p2 ) == CGAL::EQUAL );
 assert( CGAL::compare_x( p1, p4 ) == CGAL::SMALLER );
 assert( CGAL::compare_x( p4, p1 ) == CGAL::LARGER );
 assert( CGAL::compare_x( p6, p5 ) == CGAL::LARGER );

 assert( CGAL::compare_y( p3, p6 ) == CGAL::EQUAL );
 assert( CGAL::compare_y( p5, p7 ) == CGAL::SMALLER );
 assert( CGAL::compare_y( p3, p4 ) == CGAL::SMALLER );
 assert( CGAL::compare_y( p2, p1 ) == CGAL::LARGER );

 std::cout <<'.';

 CGAL::Vector_2<R>     v( RT(10), RT(20), RT(-5) );
 CGAL::Direction_2<R>  dir1(RT(-11),RT( 13)),
                      dir2(RT( 14),RT(-22)),
                      dir3(RT(-12),RT( 47)),
                      dir4(RT(- 7),RT(- 8)),
                      dir5(RT( -3), RT( 4));

 CGAL::Aff_transformation_2<R> rotate1(CGAL::ROTATION,dir1, RT(1), RT(100)),
                              rotate2(CGAL::ROTATION,dir2, RT(1), RT(100)),
                              rotate3(CGAL::ROTATION,dir3, RT(1), RT(100)),
                              rotate4(CGAL::ROTATION,dir4, RT(1), RT(100)),
                              rotate5(CGAL::ROTATION,dir5, RT(1), RT(100));

 CGAL::Point_2<R> p0( RT(1), RT(0) );

 p1 = p0.transform(rotate1);
 p2 = p0.transform(rotate2);
 p3 = p0.transform(rotate3);
 p4 = p0.transform(rotate4);
 p5 = p0.transform(rotate5);

 using CGAL::testsuite::approx_equal;
 using CGAL::testsuite::Direction_2_tag;

 const bool nonexact = boost::is_same<FT, double>::value;

 assert( approx_equal((p5 - CGAL::ORIGIN).direction(), dir5, Direction_2_tag()) );

 assert( CGAL::side_of_bounded_circle(p1, p2, p3, CGAL::Point_2<R>(CGAL::ORIGIN))\
                                      == CGAL::ON_BOUNDED_SIDE );
 assert( CGAL::side_of_bounded_circle(p1+v, p2+v, p3+v, CGAL::ORIGIN + v) \
                                      == CGAL::ON_BOUNDED_SIDE );
 assert( CGAL::side_of_bounded_circle(p1+v, p2+v, p3+v, CGAL::ORIGIN - v) \
                                      == CGAL::ON_UNBOUNDED_SIDE || nonexact);
 assert( CGAL::side_of_bounded_circle(p1, p2, p3, p4) \
                                      == CGAL::ON_BOUNDARY || nonexact);
 assert( CGAL::side_of_bounded_circle(p1+v, p2+v, p3+v, p4+v) \
                                      == CGAL::ON_BOUNDARY || nonexact);
 assert( CGAL::side_of_bounded_circle(p1+v, p3+v, p4+v, p2+v) \
                                      == CGAL::ON_BOUNDARY || nonexact);
 assert( CGAL::side_of_bounded_circle(p2+v, p4+v, p1+v, p3+v) \
                                      == CGAL::ON_BOUNDARY || nonexact);

 assert( CGAL::orientation( p1, p2, p3 ) == CGAL::POSITIVE );

 assert( CGAL::side_of_oriented_circle(p1,p2,p3,CGAL::Point_2<R>(CGAL::ORIGIN))\
                                      == CGAL::ON_POSITIVE_SIDE );
 assert( CGAL::side_of_oriented_circle(p1+v, p2+v, p3+v, CGAL::ORIGIN + v) \
                                      == CGAL::ON_POSITIVE_SIDE );
 assert( CGAL::side_of_oriented_circle(p1+v, p3+v, p2+v, CGAL::ORIGIN + v) \
                                      == CGAL::ON_NEGATIVE_SIDE );
 assert( CGAL::side_of_oriented_circle(p1+v, p2+v, p3+v, CGAL::ORIGIN - v) \
                                      == CGAL::ON_NEGATIVE_SIDE );
 assert( CGAL::side_of_oriented_circle(p2+v, p1+v, p3+v, CGAL::ORIGIN - v) \
                                      == CGAL::ON_POSITIVE_SIDE );
 assert( CGAL::side_of_oriented_circle(p1, p2, p3, p4) \
                                      == CGAL::ON_ORIENTED_BOUNDARY || nonexact);
 assert( CGAL::side_of_oriented_circle(p1+v, p2+v, p3+v, p4+v) \
                                      == CGAL::ON_ORIENTED_BOUNDARY || nonexact);
 assert( CGAL::side_of_oriented_circle(p1+v, p3+v, p4+v, p2+v) \
                                      == CGAL::ON_ORIENTED_BOUNDARY || nonexact);

 CGAL::Point_2<R> p10( RT(100), RT(100), RT(10) );
 CGAL::Point_2<R> p11( RT(-100), RT(-100), RT(10) );
 CGAL::Point_2<R> pt1( RT(-100), RT(100), RT(10) );
 CGAL::Point_2<R> pt2( CGAL::ORIGIN );
 CGAL::Point_2<R> pt3( RT(1000), RT(1000), RT(1) );

 assert( CGAL::side_of_bounded_circle(p10, p11, pt1) == CGAL::ON_BOUNDARY);
 assert( CGAL::side_of_bounded_circle(p10, p11, pt2) == CGAL::ON_BOUNDED_SIDE);
 assert( CGAL::side_of_bounded_circle(p10, p11, pt3) == CGAL::ON_UNBOUNDED_SIDE);

 // Now test squared_radius().

 assert( CGAL::squared_radius(p10, p11, pt1) == FT(200));
 assert( CGAL::squared_radius(p10, p11) == FT(200));
 assert( CGAL::squared_radius(p10) == FT(0));
 std::cout << '.';

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FURTHER_FCT_POINT_2_H
