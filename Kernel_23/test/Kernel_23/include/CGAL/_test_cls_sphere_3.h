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


#ifndef CGAL__TEST_CLS_SPHERE_3_H
#define CGAL__TEST_CLS_SPHERE_3_H

#include <CGAL/Bbox_3.h>
#include <cassert>

template <class R>
bool
_test_cls_sphere_3(const R& )
{
 std::cout << "Testing class Sphere_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;
 typename R::Sphere_3  ic;
 CGAL::Sphere_3<R> c0;

 const bool nonexact = std::is_floating_point<FT>::value;
 RT n0 =  0;
 RT n1 = 16;
 RT n2 = -4;
 RT n3 =  2;
 RT n4 =  5;
 RT n5 = 10;
 RT n6 =  6;
 RT n8 =  9;
 RT n9 = 36;
 RT n10=  3;
 RT n11=-11;
 RT n12=200;
 RT n13= 13;

 CGAL::Point_3<R> p0( n1, n2, n1, -n2);  // ( 4, -1, 4)
 CGAL::Point_3<R> p1( n6, n8, n8, n10);  // ( 2,  3, 3)
 CGAL::Point_3<R> p2( n2, n0, n0,  n2);  // ( 1,  0, 0)
 CGAL::Point_3<R> p3( n5, n5, n5,  n4);  // ( 2,  2, 2)
 CGAL::Point_3<R> p4( n0, n0, n2,  n2);  // ( 0,  0, 1)
 CGAL::Point_3<R> p5( n8, n6, n10, n10); // ( 3,  2, 1)

 CGAL::Vector_3<R> vx = p2 - CGAL::ORIGIN;
 CGAL::Vector_3<R> vy = p4 - CGAL::ORIGIN;
 CGAL::Vector_3<R> v1 = p1 - CGAL::ORIGIN;

 CGAL::Sphere_3<R> c1( p0, p1, p2, p4);
 CGAL::Sphere_3<R> c2( p0, p1, p3, p5);
 CGAL::Sphere_3<R> c3( p1, p0, p2, p4);
 const CGAL::Sphere_3<R> c4( p3, FT( n9 ));    // n9 = (n6)^2
 CGAL::Vector_3<R> vx6 = vx * n6;
 CGAL::Vector_3<R> vy6 = vy * n6;
 CGAL::Sphere_3<R> c5( p3 - vx6, p3 + vx6, p3 + vy6);
 CGAL::Sphere_3<R> c6( c3 );
 CGAL::Sphere_3<R> c7( p3, n9, CGAL::POSITIVE);
 CGAL::Sphere_3<R> c8( p3, n9, CGAL::NEGATIVE);
 CGAL::Sphere_3<R> cc( p3 - vx6, p3 + vx6);
 CGAL::Sphere_3<R> cp( p3 - vx6, p3 + vx6, CGAL::POSITIVE);
 CGAL::Sphere_3<R> cn( p3 - vx6, p3 + vx6, CGAL::NEGATIVE);
 CGAL::Sphere_3<R> cc3( p3 - vx6, p3 + vx6, p3 + vy6);
 CGAL::Sphere_3<R> cp3( p3 - vx6, p3 + vx6, p3 + vy6, CGAL::POSITIVE);
 CGAL::Sphere_3<R> cn3( p3 - vx6, p3 + vx6, p3 + vy6, CGAL::NEGATIVE);
 c0 = c3;

 assert( c1 == c1 );
 assert( c1 != c2 );
 assert( c3 == c0 );
 assert( c0 == c3 );
 assert( c3 == c6 );
 assert( c7 != c8 );
 assert( c4 == c7 );
 assert( cc == cp );
 assert( cn != cp );
 assert( cc != c8 );
 assert( cc == c7 );

 assert( c5.center() == p3 || nonexact);
 assert( cc.center() == p3 );
 assert( c5.squared_radius() == FT( n9 ) );
 assert( c4.squared_radius() == cc.squared_radius() );
 assert( c4 == c5 || nonexact);
 assert( c4 == c7 );
 assert( c4 != c8 );
 assert( cn == cp.opposite() );
 assert( cn3 == cp3.opposite() );
 assert( c7.opposite() == c8 );
 assert( c8.opposite() == c7 );
 assert( c1.opposite() == c3 || nonexact );
 assert( c3.opposite() == c1 || nonexact );
 assert( c7.orientation() == CGAL::POSITIVE );
 assert( c8.orientation() == CGAL::NEGATIVE );
 assert( c5.orientation() == CGAL::POSITIVE );
 assert( cc.orientation() == CGAL::POSITIVE );
 assert( cp.orientation() == CGAL::POSITIVE );
 assert( cn.orientation() == CGAL::NEGATIVE );
 assert( cc3.orientation() == CGAL::POSITIVE );
 assert( cp3.orientation() == CGAL::POSITIVE );
 assert( cn3.orientation() == CGAL::NEGATIVE );

 std::cout << '.';

 assert( c4.center() == p3 );
 assert( c5.center() == p3 || nonexact);
 assert( c4.squared_radius() == FT( n9 ) );
 assert( c5.squared_radius() == FT( n9 ) );
 assert( c8.squared_radius() == FT( n9 ) );

 assert( c7.bounded_side( p3 + vx*n2 ) == CGAL::ON_BOUNDED_SIDE );
 assert( c7.bounded_side( p3 + vy*n11 ) == CGAL::ON_UNBOUNDED_SIDE );
 assert( c7.bounded_side( p3 - vy6 ) == CGAL::ON_BOUNDARY );
 assert( c8.bounded_side( p3 + vx*n2 ) == CGAL::ON_BOUNDED_SIDE );
 assert( c8.bounded_side( p3 + vy*n11 ) == CGAL::ON_UNBOUNDED_SIDE );
 assert( c8.bounded_side( p3 - vy6 ) == CGAL::ON_BOUNDARY );
 assert( cc.has_on_boundary( p3 + vy6) );
 assert( cc.has_on_boundary( p3 - vx6) );

 std::cout << '.';

 RT sin1;
 RT cos1;
 RT den1;
 RT sin2;
 RT cos2;
 RT den2;
 RT sin3;
 RT cos3;
 RT den3;
 RT sin4;
 RT cos4;
 RT den4;
 RT sin5;
 RT cos5;
 RT den5;
 RT RT0 = RT(0);
 CGAL::rational_rotation_approximation<RT>(n11,n13, sin1, cos1, den1, -n2,n12);
 CGAL::rational_rotation_approximation<RT>(-n8, n9, sin2, cos2, den2, -n2,n12);
 CGAL::rational_rotation_approximation<RT>( n5,-n1, sin3, cos3, den3, -n2,n12);
 CGAL::rational_rotation_approximation<RT>(-n5,-n11,sin4, cos4, den4, -n2,n12);
 CGAL::rational_rotation_approximation<RT>(-n2, n2, sin5, cos5, den5, -n2,n12);

 CGAL::Aff_transformation_3<R> rotate1( sin1, cos1,  RT0,  RT0,
                                       -cos1, sin1,  RT0,  RT0,
                                         RT0,  RT0, den1,  RT0,
                                                          den1 );

 CGAL::Aff_transformation_3<R> rotate2( sin2, cos2,  RT0,  RT0,
                                       -cos2, sin2,  RT0,  RT0,
                                         RT0,  RT0, den2,  RT0,
                                                          den2 );

 CGAL::Aff_transformation_3<R> rotate3( sin3,  RT0, cos3,  RT0,
                                         RT0, den3,  RT0,  RT0,
                                       -cos3,  RT0, sin3,  RT0,
                                                          den3 );

 CGAL::Aff_transformation_3<R> rotate4( den4,  RT0,  RT0,  RT0,
                                         RT0, sin4, cos4,  RT0,
                                         RT0,-cos4, sin4,  RT0,
                                                          den4 );

 CGAL::Aff_transformation_3<R> rotate5 = rotate1 * rotate4;

 CGAL::Point_3<R> ori = CGAL::Point_3<R>( RT(0), RT(0), RT(0));
 CGAL::Point_3<R> p6  = p2.transform( rotate1 );
 assert( CGAL::compare_distance_to_point( ori, p2, p6) == CGAL::EQUAL || nonexact );
 CGAL::Point_3<R> p7  = p2.transform( rotate2 );
 assert( CGAL::compare_distance_to_point( ori, p2, p7) == CGAL::EQUAL || nonexact );
 CGAL::Point_3<R> p8  = p2.transform( rotate3 );
 assert( CGAL::compare_distance_to_point( ori, p2, p8) == CGAL::EQUAL || nonexact );
 CGAL::Point_3<R> p9  = p2.transform( rotate4 );
 assert( CGAL::compare_distance_to_point( ori, p2, p9) == CGAL::EQUAL );
 CGAL::Point_3<R> p10 = p2.transform( rotate5 );
 assert( CGAL::compare_distance_to_point( ori, p2, p10) == CGAL::EQUAL || nonexact );
 p6 = p6 + v1;
 p7 = p7 + v1;
 p8 = p8 + v1;
 p9 = p9 + v1;
 p10 = p10 + v1;
 CGAL::Sphere_3<R> c10 (p6, p8, p7, p9);
 assert( c10.center() == ori + v1 || nonexact );
 assert( c10.orientation() == CGAL::POSITIVE );
 assert( c10.opposite().orientation() == CGAL::NEGATIVE );

 assert( c10.oriented_side(c10.center() ) == CGAL::ON_POSITIVE_SIDE );
 assert( c10.oriented_side(CGAL::ORIGIN + v1 + vx/n2 ) \
         == CGAL::ON_POSITIVE_SIDE );
 assert( c10.oriented_side(CGAL::ORIGIN + v1 + vx*n2 ) \
         == CGAL::ON_NEGATIVE_SIDE );
 assert( c10.oriented_side(p9 ) == CGAL::ON_ORIENTED_BOUNDARY || nonexact );
 assert( c10.has_on_boundary(p9) || nonexact );
 assert( c10.has_on_boundary(p4 + v1) || nonexact );
 CGAL::Point_3<R> p11( n4, n4, n4, n3) ; // (2.5, 2.5, 2.5)
 CGAL::Point_3<R> p12( n5, n5, n5, n3) ; // ( 5 ,  5,   5 )
 assert( c10.has_on_bounded_side( p11 ) );
 assert( ! c10.has_on_bounded_side( p12 ) );
 assert( c10.has_on_unbounded_side( p12 ) );
 assert( c10.has_on_positive_side( p11 ) );
 assert( c10.has_on_negative_side( p12 ) );
 assert( c10.opposite().has_on_negative_side( p11 ) );
 assert( c10.opposite().has_on_positive_side( p12 ) );
 assert( c10.has_on_boundary( p6 ) || nonexact );
 assert( c10.has_on_boundary( p8 ) || nonexact );


 std::cout << '.';

 CGAL::Sphere_3<R> c11( p0 );
 CGAL::Sphere_3<R> c12( p0, CGAL::POSITIVE );
 CGAL::Sphere_3<R> c13( p0, CGAL::NEGATIVE );
 assert( c11.orientation() == CGAL::POSITIVE );
 assert( c12.orientation() == CGAL::POSITIVE );
 assert( c13.orientation() == CGAL::NEGATIVE );
 assert( c11.is_degenerate() );
 assert( c12.is_degenerate() );
 assert( c13.is_degenerate() );

 std::cout << '.';

 CGAL::Bbox_3 bb = c4.bbox();
 assert(bb.xmin() <= -4.0);
 assert(bb.xmax() >= 8.0);
 assert(bb.ymin() <= -4.0);
 assert(bb.ymax() >= 8.0);
 assert(bb.zmin() <= -4.0);
 assert(bb.zmax() >= 8.0);

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_SPHERE_3_H
