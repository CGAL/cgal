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
// file          : _test_cls_ray_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_CLS_RAY_2_C
#define CGAL__TEST_CLS_RAY_2_C
#ifndef CGAL__TEST_CLS_RAY_2_H
#include <CGAL/_test_cls_ray_2.h>
#endif // CGAL__TEST_CLS_RAY_2_H

template <class R>
bool
_test_cls_ray_2(const R& )
{
 std::cout << "Testing class Ray_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Ray_2  ir;
 CGAL::Ray_2<R> r0;

 RT  n0 = 0;
 RT  n2 = 2;
 RT  n3 = 3;
 RT  n4 = 4;
 RT  n5 = 5;
 RT  n8 = 8;
 RT  n9 = 9;
 RT n10 =10;

 CGAL::Point_2<R> p1( n2, n8, n2);
 CGAL::Point_2<R> p2( n10, n4, n2);
 CGAL::Point_2<R> p3( n9, n9, n3);
 CGAL::Point_2<R> p4( n10, n8, n2);
 CGAL::Direction_2<R> d12( p2 - p1);
 CGAL::Direction_2<R> d24( p4 - p2);

 CGAL::Ray_2<R> r1( p1, p2);
 CGAL::Ray_2<R> r2( p1, d12);
 CGAL::Ray_2<R> r3( p3, p2);
 CGAL::Ray_2<R> r4( p2, d24);
 CGAL::Ray_2<R> r5( p2, p4);
 CGAL::Ray_2<R> r6( p1, p4);
 r0 = r3;

 std::cout << '.';

 assert( r1 == r1 );
 assert( r2 == r1 );
 assert( r4 == r5 );
 assert( r0 == r3 );
 assert( r6 != r1 );
 assert( r1 != r3 );
 assert( r1 != r5 );

 assert( r2.source() == p1 );
 assert( r0.source() == p3 );
 assert( r4.source() == r4.point(0) );

 assert( r1.direction() == d12 );
 assert( r2.direction() == d12 );
 assert( r3.direction() == r1.direction() );

 std::cout << '.';

 assert( r3.supporting_line() == r1.supporting_line() );
 assert( r5.supporting_line() == CGAL::Line_2<R>( p2, p4 ) );

 assert( r4.opposite() == CGAL::Ray_2<R>( p2, -d24 ) );
 assert( r1.opposite() == CGAL::Ray_2<R>( p1, -d12 ) );
 assert( r2.opposite().opposite() == r2 );

 assert( r6.is_horizontal() );
 assert( ! r0.is_horizontal() );
 assert( r5.is_vertical() );
 assert( ! r5.is_horizontal() );

 assert( r1.has_on( p1 ) );
 assert( r1.has_on( p2 ) );
 assert( r1.has_on( p3 ) );
 assert( r3.opposite().has_on( p1 ) );
 assert( ! r1.has_on( p4 ) );
 assert( ! r0.has_on( CGAL::Point_2<R>( n8, n5, n8 )) );
 assert( r4.has_on( r4.point(7)) );
 assert( r3.collinear_has_on( r3.point(7)) );
 assert( r1.collinear_has_on( p3) );
 assert( ! r3.collinear_has_on( p1 ) );

 std::cout << '.';

 assert( CGAL::Ray_2<R>( p1, p1).is_degenerate() );
 assert( ! r0.is_degenerate() );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_RAY_2_C
