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
// file          : _test_cls_direction_3.C
// revision      : 3.8
// revision_date : 08 Oct 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_CLS_DIRECTION_3_C
#define CGAL__TEST_CLS_DIRECTION_3_C
#ifndef CGAL__TEST_CLS_DIRECTION_3_H
#include <CGAL/_test_cls_direction_3.h>
#endif // CGAL__TEST_CLS_DIRECTION_3_H

template <class R>
bool
_test_cls_direction_3(const R& )
{
 std::cout << "Testing class Direction_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Direction_3  id;

 CGAL::Direction_3<R> d0;
 CGAL::Direction_3<R> d1(id);

 std::cout << '.';
 RT   n0 = 10;
 RT  n1 = 8;
 RT  n2 = 4;
 RT  n3 = 2;

 CGAL::Vector_3<R>  v( n1, n2, n3);
 CGAL::Direction_3<R> d2(v);
 CGAL::Direction_3<R> d3( n0, n1, n2);
 CGAL::Direction_3<R> d4( d3 );
 CGAL::Direction_3<R> d5 = d3;

 assert( d3 == d3 );
 assert( d3 == d4 );
 assert( d5 == d3 );
 assert( d2 != d3 );
 assert( d3 != d2 );

 std::cout << '.';
#if (__GNUG__ == 2) && (__GNUC_MINOR__==91)
 CGAL::Vector_3<R> vv = d2.to_vector();
#else
 CGAL::Vector_3<R> vv = d2.vector();
#endif // egcs 2.91.66
 assert( v == vv );

 d0 = -d3;

 assert( d0 != d3 );
 assert( d3 == -d0);

 std::cout << '.';
 assert( d3.delta(0) == n0 );
 assert( d3.delta(1) == n1 );
 assert( d3.delta(2) == n2 );
 assert( d3.delta(0) == d3.dx() );
 assert( d3.delta(1) == d3.dy() );
 assert( d3.delta(2) == d3.dz() );

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_DIRECTION_3_C
