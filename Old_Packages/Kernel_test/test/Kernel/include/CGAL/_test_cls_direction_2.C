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
// file          : _test_cls_direction_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_CLS_DIRECTION_2_C
#define CGAL__TEST_CLS_DIRECTION_2_C
#ifndef CGAL__TEST_CLS_DIRECTION_2_H
#include <CGAL/_test_cls_direction_2.h>
#endif // CGAL__TEST_CLS_DIRECTION_2_H

template <class R>
bool
_test_cls_direction_2(const R& )
{
 std::cout << "Testing class Direction_2" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Direction_2  id;
 CGAL::Direction_2<R> d0;
 CGAL::Direction_2<R> d1(id);

 std::cout << '.';
 RT  n0 = 10;
 RT  n1 = 8;
 RT  n2 = 4;
 RT  n3 = 2;

 CGAL::Vector_2<R>  v( n1, n2);      // (8,4)
 CGAL::Direction_2<R> d2(v);
 CGAL::Direction_2<R> d3( n0, n1);   // (10,8)
 CGAL::Direction_2<R> d4( d3 );
 CGAL::Direction_2<R> d5 = d3;
 CGAL::Point_2<R> p1(n3, n2); // (4,2)
 CGAL::Line_2<R> l1(p1, d3);
 CGAL::Ray_2<R> r1(p1, d2);
 CGAL::Segment_2<R> s1(p1, p1 + d4.vector());
 CGAL::Direction_2<R> d6( l1);
 CGAL::Direction_2<R> d7( r1);
 CGAL::Direction_2<R> d8( s1);

 assert( d3 == d3 );
 assert( d3 == d4 );
 assert( d5 == d3 );
 assert( d2 != d3 );
 assert( d3 != d2 );
 assert( d6 == d8 );
 assert( d6 != d7 );
 assert( d7 != d8 );

 std::cout << '.';
 CGAL::Vector_2<R> vv = d2.vector();
 assert( v == vv );

 d0 = -d3;

 assert( d0 != d3 );
 assert( d3 == -d0);

 std::cout << '.';
 assert( d3.delta(0) == n0 );
 assert( d3.delta(1) == n1 );
 assert( d3.delta(0) == d3.dx() );
 assert( d3.delta(1) == d3.dy() );

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_DIRECTION_2_C
