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
// file          : _test_fct_point_vector_2.C
// revision      : 2.1
// revision_date : 05 Aug 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL__TEST_FCT_POINT_VECTOR_2_C
#define CGAL__TEST_FCT_POINT_VECTOR_2_C
#ifndef CGAL__TEST_FCT_POINT_VECTOR_2_H
#include <CGAL/_test_fct_point_vector_2.h>
#endif // CGAL__TEST_FCT_POINT_VECTOR_2_H

template <class R>
bool
_test_fct_point_vector_2(const R& )
{
 std::cout << "Testing functions Point_2 Vector_2" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 RT  n1( 12 );
 RT  n2( -4 );
 RT  n3(  6 );
 RT  n4(  2 );
 RT  n5(  9 );
 RT  n6(-18 );
 RT  n7( 36 );
 RT  n8(  3 );
 RT  n9( 15 );
 RT n10( -8 );

 CGAL::Vector_2<R>  v0(CGAL::NULL_VECTOR);
 CGAL::Vector_2<R>  v1(n1, n2, n4);
 CGAL::Vector_2<R>  v2(n5, n6, n8);
 CGAL::Vector_2<R>  v3(n5, n10);
 CGAL::Vector_2<R>  v4(n8, -n2);

 std::cout << '.';

 CGAL::Point_2<R> p0(CGAL::ORIGIN);
 CGAL::Point_2<R> p1 = CGAL::ORIGIN + v1;
 CGAL::Point_2<R> p2 = CGAL::ORIGIN + v2;
 CGAL::Point_2<R> p3 = CGAL::ORIGIN + v3;

 assert( CGAL::ORIGIN + v2 == CGAL::Point_2<R>( n5, n6, n8) );
 assert( CGAL::ORIGIN - v2 == CGAL::Point_2<R>( -n5, -n6, n8) );
 assert( p1 - p1 == v0 );
 assert( p1 - p0 == p1 - CGAL::ORIGIN);
 assert( p1 - p2 == v4 );
 assert( p2 + v4 == p1 );
 assert( p3 - v1 == p2 );
 assert( p3 - p1 == v2 );

 std::cout << "..";
 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_POINT_VECTOR_2_C
