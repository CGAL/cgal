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
// file          : _test_cls_triangle_3.C
// revision      : 2.0.5
// revision_date : 24 Mar 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL__TEST_CLS_TRIANGLE_3_C
#define CGAL__TEST_CLS_TRIANGLE_3_C
#ifndef CGAL__TEST_CLS_TRIANGLE_3_H
#include <CGAL/_test_cls_triangle_3.h>
#endif // CGAL__TEST_CLS_TRIANGLE_3_H

template <class R>
bool
_test_cls_triangle_3(const R& )
{
 cout << "Testing class Triangle_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Triangle_3 it;
 CGAL::Triangle_3<R> t0(it);

 RT n0 =  0;
 RT n1 = 12;
 RT n2 = 16;
 RT n3 = -4;
 RT n4 =  2;
 RT n5 =  3;
 RT n6 = 30;
 RT n7 =  9;
 RT n8 = 24;
 RT n9 =  8;

 CGAL::Point_3<R> p1( n1, n2, n3, n4);
 CGAL::Point_3<R> p2( n2, n9, n3,-n3);
 CGAL::Point_3<R> p3( n5, n6, n1, n5);
 CGAL::Point_3<R> p4( n7, n7, n8, n5);

 CGAL::Point_3<R> ps3( n0, n0, n7, n5);
 CGAL::Point_3<R> ps2( n0, n7, n0, n5);
 CGAL::Point_3<R> ps1( n7, n0, n0, n5);

 CGAL::Triangle_3<R> t1(p1,p2,p3);
 CGAL::Triangle_3<R> t2(p4,p2,p3);
 CGAL::Triangle_3<R> t3(ps1,ps2,ps3);
 CGAL::Triangle_3<R> t4(ps2,ps1,ps3);
 CGAL::Triangle_3<R> t5( t1 );
 t0 = t1;

 assert( t0 == t0 );
 assert( t0 == t1 );
 assert( t5 == t1 );
 assert( t2 != t4 );
 assert( t3 != t4 );

 cout <<'.';

 CGAL::Plane_3<R> pl1( p1,p2,p3);
 CGAL::Plane_3<R> pl2( p4,p2,p3);
 assert( t1.supporting_plane() == pl1 );
 assert( t2.supporting_plane() == pl2 );
 assert( t3.supporting_plane() == t4.supporting_plane().opposite() );

 cout <<'.';

 assert( t1.has_on(p3) );
 assert( t1.has_on(p2) );
 assert( t2.has_on(p4) );
 assert( ! t1.has_on(p4) );
 CGAL::Point_3<R> pt( n7, n7, n7, n7);
 assert( t3.has_on( pt ) );
 assert( t4.has_on( pt ) );

 assert( t1.vertex(0) == p1 );
 assert( t1.vertex(1) == p2 );
 assert( t1.vertex(2) == p3 );
 assert( t4[0] == ps2 );
 assert( t4[1] == ps1 );
 assert( t4[2] == ps3 );

 cout <<'.';

 CGAL::Triangle_3<R> tdeg1( p3,p3,p1);
 CGAL::Triangle_3<R> tdeg2( p3,p3,p3);
 assert( tdeg1.is_degenerate() );
 assert( tdeg2.is_degenerate() );

 cout << "done" << endl;
 return true;
}

#endif // CGAL__TEST_CLS_TRIANGLE_3_C
