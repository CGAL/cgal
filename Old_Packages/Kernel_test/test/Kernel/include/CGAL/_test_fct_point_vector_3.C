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
// file          : _test_fct_point_vector_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_FCT_POINT_VECTOR_3_C
#define CGAL__TEST_FCT_POINT_VECTOR_3_C

#include <CGAL/_test_fct_point_vector_3.h>

template <class R>
bool
_test_fct_point_vector_3(const R& )
{
 std::cout << "Testing functions Point_3 Vector_3" ;

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

 CGAL::Vector_3<R>  v0(CGAL::NULL_VECTOR);
 CGAL::Vector_3<R>  v1(n1, n2, n3, n4);
 CGAL::Vector_3<R>  v2(n5, n6, n7, n8);
 CGAL::Vector_3<R>  v3(n5, n10, n9);
 CGAL::Vector_3<R>  v4(n8, -n2, -n5);

 std::cout << '.';

 CGAL::Point_3<R> p0(CGAL::ORIGIN);
 CGAL::Point_3<R> p1 = CGAL::ORIGIN + v1;
 CGAL::Point_3<R> p2 = CGAL::ORIGIN + v2;
 CGAL::Point_3<R> p3 = CGAL::ORIGIN + v3;

 CGAL::Vector_3<R>  v5(p1, p2);

 assert( CGAL::ORIGIN + v2 == CGAL::Point_3<R>( n5, n6, n7, n8) );
 assert( CGAL::ORIGIN - v2 == CGAL::Point_3<R>( -n5, -n6, -n7, n8) );
 assert( p1 - p1 == v0 );
 assert( p1 - p0 == p1 - CGAL::ORIGIN);
 assert( p1 - p2 == v4 );
 assert( p2 + v4 == p1 );
 assert( p3 - v1 == p2 );
 assert( p3 - p1 == v2 );
 assert( p2 - p1 == v5 );

 std::cout << "..";
 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_POINT_VECTOR_3_C
