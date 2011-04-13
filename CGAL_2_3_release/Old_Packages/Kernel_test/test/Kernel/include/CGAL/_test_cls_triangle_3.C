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
// file          : _test_cls_triangle_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
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
 std::cout << "Testing class Triangle_3" ;

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

 CGAL::Point_3<R> p1( n1, n2, n3, n4);  // (6, 8, -2)
 CGAL::Point_3<R> p2( n2, n9, n3,-n3);  // (-4, 2, -1)
 CGAL::Point_3<R> p3( n5, n6, n1, n5);  // (1, 10, 4)
 CGAL::Point_3<R> p4( n7, n7, n8, n5);  // (3, 3, 6)
 CGAL::Point_3<R> p5( n2, n1, n0, n4);  // (8, 6, 0)
 CGAL::Point_3<R> p6( n0, n0, n1, n5);  // (0, 0, 4)

 CGAL::Point_3<R> ps3( n0, n0, n7, n5); // (0, 0, 3)
 CGAL::Point_3<R> ps2( n0, n7, n0, n5); // (0, 3, 0)
 CGAL::Point_3<R> ps1( n7, n0, n0, n5); // (3, 0, 0)

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

 std::cout <<'.';

 CGAL::Plane_3<R> pl1( p1,p2,p3);
 CGAL::Plane_3<R> pl2( p4,p2,p3);
 assert( t1.supporting_plane() == pl1 );
 assert( t2.supporting_plane() == pl2 );
 assert( t3.supporting_plane() == t4.supporting_plane().opposite() );

 std::cout <<'.';

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

 std::cout <<'.';

 CGAL::Triangle_3<R> tdeg1( p3,p3,p1);
 CGAL::Triangle_3<R> tdeg2( p3,p3,p3);
 assert( tdeg1.is_degenerate() );
 assert( tdeg2.is_degenerate() );

 std::cout <<'.';

 assert ( tdeg1.squared_area() == FT(0) );
 assert ( tdeg2.squared_area() == FT(0) );
 assert ( t5.squared_area() == t1.squared_area() );
 assert ( t3.squared_area() == t4.squared_area() );
 CGAL::Triangle_3<R> t6(ps3,p5,p6);
 assert ( t6.squared_area() == FT(25) );

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_TRIANGLE_3_C
