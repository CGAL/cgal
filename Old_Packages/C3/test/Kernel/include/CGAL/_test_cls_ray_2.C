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
// source        : test_kernel_2.fw
// file          : _test_cls_ray_2.C
// revision      : 2.0.5
// revision_date : 24 Mar 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL__TEST_CLS_RAY_2_C
#define CGAL__TEST_CLS_RAY_2_C
#ifndef CGAL__TEST_CLS_RAY_2_H
#include <CGAL/_test_cls_ray_2.h>
#endif // CGAL__TEST_CLS_RAY_2_H

template <class gnuR>
bool
_test_cls_ray_2(const gnuR& )
{
 cout << "Testing class Ray_2";

 typedef typename  gnuR::RT    RT;
 typedef typename  gnuR::FT    FT;

 typename gnuR::Ray_2  ir;
 CGAL::Ray_2<gnuR> r0;

 RT  n0 = 0;
 RT  n2 = 2;
 RT  n3 = 3;
 RT  n4 = 4;
 RT  n5 = 5;
 RT  n8 = 8;
 RT  n9 = 9;
 RT n10 =10;

 CGAL::Point_2<gnuR> p1( n2, n8, n2);
 CGAL::Point_2<gnuR> p2( n10, n4, n2);
 CGAL::Point_2<gnuR> p3( n9, n9, n3);
 CGAL::Point_2<gnuR> p4( n10, n8, n2);
 CGAL::Direction_2<gnuR> d12( p2 - p1);
 CGAL::Direction_2<gnuR> d24( p4 - p2);

 CGAL::Ray_2<gnuR> r1( p1, p2);
 CGAL::Ray_2<gnuR> r2( p1, d12);
 CGAL::Ray_2<gnuR> r3( p3, p2);
 CGAL::Ray_2<gnuR> r4( p2, d24);
 CGAL::Ray_2<gnuR> r5( p2, p4);
 CGAL::Ray_2<gnuR> r6( p1, p4);
 r0 = r3;

 cout << '.';

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

 cout << '.';

 assert( r3.supporting_line() == r1.supporting_line() );
 assert( r5.supporting_line() == CGAL::Line_2<gnuR>( p2, p4 ) );

 assert( r4.opposite() == CGAL::Ray_2<gnuR>( p2, -d24 ) );
 assert( r1.opposite() == CGAL::Ray_2<gnuR>( p1, -d12 ) );
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
 assert( ! r0.has_on( CGAL::Point_2<gnuR>( n8, n5, n8 )) );
 assert( r4.has_on( r4.point(7)) );
 assert( r3.collinear_has_on( r3.point(7)) );
 assert( r1.collinear_has_on( p3) );
 assert( ! r3.collinear_has_on( p1 ) );

 cout << '.';

 assert( CGAL::Ray_2<gnuR>( p1, p1).is_degenerate() );
 assert( ! r0.is_degenerate() );

 cout << "done" << endl;
 return true;
}
#endif // CGAL__TEST_CLS_RAY_2_C
