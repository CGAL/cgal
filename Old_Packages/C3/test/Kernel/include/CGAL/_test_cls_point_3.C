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
// source        : test_kernel_3.fw
// file          : _test_cls_point_3.C
// revision      : 2.0.5
// revision_date : 24 Mar 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL__TEST_CLS_POINT_3_C
#define CGAL__TEST_CLS_POINT_3_C
#ifndef CGAL__TEST_CLS_POINT_3_H
#include <CGAL/_test_cls_point_3.h>
#endif // CGAL__TEST_CLS_POINT_3_H

template <class R>
bool
_test_cls_point_3(const R& )
{
 cout << "Testing class Point_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Point_3       ip;

 CGAL::Point_3<R>  p1;
 CGAL::Point_3<R>  p2(ip);
 CGAL::Point_3<R>  p0(CGAL::ORIGIN);

 RT  n1(-35 );
 RT  n2( 50 );
 RT  n3(-20 );
 RT  n4(  5 );

 CGAL::Point_3<R>  p3(n1, n2, n3);
 CGAL::Point_3<R>  p4(n1, n2, n3, n4);
 CGAL::Point_3<R>  p5(n1, n2, n3, n4);
 CGAL::Point_3<R>  p6( p5 );
                  p1 = p4;

 cout << '.';

 assert( p4 == p5 );
 assert( p5 == p6 );
 assert( p4 == p6 );
 assert( p1 == p6 );

 assert( p3 != p4 );
 assert( p0 != p1 );

 assert( p0 == CGAL::ORIGIN);
 assert( p1 != CGAL::ORIGIN);
 // Doesn't work; Point_2::operator== can't be used :(
#ifdef ENHANCED
 assert( CGAL::ORIGIN == p0 );
 assert( CGAL::ORIGIN != p1 );
#endif // ENHANCED

 assert( p3.hx() == n1 );   // don't replace p3
 assert( p3.hy() == n2 );
 assert( p3.hz() == n3 );

 assert( FT(p5.hx()) / FT(p5.hw()) == FT( n1) / FT( n4) );
 assert( FT(p5.hy()) / FT(p5.hw()) == FT( n2) / FT( n4) );
 assert( FT(p5.hz()) / FT(p5.hw()) == FT( n3) / FT( n4) );

 assert( p5.x() == FT( n1) / FT( n4 ) );
 assert( p5.y() == FT( n2) / FT( n4 ) );
 assert( p5.z() == FT( n3) / FT( n4 ) );

 cout << '.';

 assert( p3.homogeneous(0) == p3.hx() );  // don't replace p3
 assert( p3.homogeneous(1) == p3.hy() );
 assert( p3.homogeneous(2) == p3.hz() );
 assert( p3.homogeneous(3) == p3.hw() );
 assert( p6.cartesian(0) == p6.x() );
 assert( p6.cartesian(1) == p6.y() );
 assert( p6.cartesian(2) == p6.z() );

 cout << '.';

 assert( p0.dimension() == 3 );

 cout << "done" << endl;
 return true;
}
#endif // CGAL__TEST_CLS_POINT_3_C
