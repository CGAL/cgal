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
// file          : _test_cls_line_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_CLS_LINE_3_C
#define CGAL__TEST_CLS_LINE_3_C
#ifndef CGAL__TEST_CLS_LINE_3_H
#include <CGAL/_test_cls_line_3.h>
#endif // CGAL__TEST_CLS_LINE_3_H

template <class R>
bool
_test_cls_line_3(const R& )
{
 std::cout << "Testing class Line_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Line_3 il;
 CGAL::Line_3<R> l0( il );
 CGAL::Line_3<R> l1;

 RT n1 =  3;
 RT n2 = 53;
 RT n3 =-14;
 RT n4 = 28;
 RT n5 = 16;
 RT n6 = -6;
 RT n7 = 11;
 RT n8 =-17;
 RT n9 = 30;

 CGAL::Point_3<R> p1( n1, n2, n3);
 CGAL::Point_3<R> p2( n4, n5, n6);
 CGAL::Point_3<R> p3( n7, n8, n9);

 CGAL::Line_3<R> l2(p1,p2);
 CGAL::Line_3<R> l3(l2);

 assert( l2 == l2);
 assert( l2 == l3);

 CGAL::Direction_3<R> dir( n9, n3, n1);
 CGAL::Line_3<R> l4(p3, dir);
 assert( l2 != l4);

 CGAL::Segment_3<R> seg(p1,p2);
 CGAL::Ray_3<R>     ray(p2,p1);
 CGAL::Line_3<R>    l5(seg);
 CGAL::Line_3<R>    l6(ray);
 assert( l2 == l5);

 std::cout <<'.';

 assert( l2 == l5 );
 assert( l2.direction() == l5.direction() );
 assert( l5.direction() ==  - l6.direction() );
 assert( l5.has_on( p1 ) );
 assert( l5.has_on( p2 ) );
 assert( l5.has_on( l5.point() ));
 assert( l6.has_on( p1 ) );
 assert( l6.has_on( p2 ) );
 assert( l6.has_on( l5.point() ));
 assert( l5.opposite() == l6 );
 assert( l2.opposite() == l6 );
 assert( l5 != l6 );

 CGAL::Plane_3<R> pl = l6.perpendicular_plane( l6.point() );
 CGAL::Plane_3<R> plstrich(l6.point(), l6.direction() );
 assert( pl == plstrich );
 assert( pl.orthogonal_direction() == l6.direction() );
 CGAL::Plane_3<R> plzweistrich(l6.point(), l5.direction() );
 assert( plzweistrich.opposite() == pl );

 std::cout << '.';

 assert( l4.point(2) - l4.point(1) == l4.point(1) - l4.point(0) );
 CGAL::Point_3<R> p1l4proj = l4.projection(p1);
 assert( l4.has_on( p1l4proj ) );
 assert( l4.perpendicular_plane( p1l4proj ).has_on( p1l4proj ) );
 assert( l4.perpendicular_plane( p1l4proj ).has_on( p1 ) );
 CGAL::Point_3<R> p4 = l4.projection(p2);
 CGAL::Point_3<R> p5 = l4.projection(p3);
 assert(  ( l4.direction() == ( p5 - p4 ).direction() )\
        ||( l4.direction() == ( p4 - p5 ).direction() )  );
 assert( l5.direction() == - l6.direction() );

 std::cout <<'.';

 assert( l2.has_on(p1) );
 assert( l2.has_on(p2) );
 assert( l4.has_on(p4) );
 assert( l4.has_on(p5) );
 assert( CGAL::Line_3<R>(p1,p1).is_degenerate() );

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_LINE_3_C
