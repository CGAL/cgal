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
// file          : _test_fct_vector_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_FCT_VECTOR_2_C
#define CGAL__TEST_FCT_VECTOR_2_C
#ifndef CGAL__TEST_FCT_VECTOR_2_H
#include <CGAL/_test_fct_vector_2.h>
#endif // CGAL__TEST_FCT_VECTOR_2_H

template <class R>
bool
_test_fct_vector_2(const R& )
{
 std::cout << "Testing functions Vector_2" ;

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
 RT n11( 24 );

 CGAL::Vector_2<R>  v0(CGAL::NULL_VECTOR);   // ( 0, 0)
 CGAL::Vector_2<R>  v1(n1, n2, n4);         // ( 6,-2)
 CGAL::Vector_2<R>  v2(n5, n6, n8);         // ( 3,-6)
 CGAL::Vector_2<R>  v3(n5, n10);            // ( 9,-8)
 CGAL::Vector_2<R>  v4(n8, -n2 );           // ( 3, 4)
 CGAL::Vector_2<R> mv4(-n8, n2 );           // (-3,-4)

 RT  n15( -5 );
 CGAL::Vector_2<R>  v8( n1,  n3, -n4);      // (-6,-3)
 CGAL::Vector_2<R>  v9( RT(0), n15);        // ( 0,-5)
 CGAL::Vector_2<R>  v10( n8,-n2, n4);       // (1.5,2)
 CGAL::Vector_2<R>  v11(-n6, n11, n8);      // ( 6, 8)


 assert( v1 + v2 == v3 );
 assert( v1 - v2 == v4 );
 assert( v3 - v1 == v2 );
 assert( v3 - v2 == v1 );
 assert( v4 + v2 == v1 );
 assert( v4 + v2 == v3 - v2);
 assert( v9 == v1 + v8 );

 std::cout << '.';

 assert( (-(- v1)) == v1 );
 assert( -v4 == mv4);
 assert( mv4 == v2 - v1);
 assert( -( v1 - v2) == mv4);
 assert( v1 + v0 == v1 );
 assert( v0 - v4 == mv4 );
 assert( v0 + v4 == v4 );
 assert( v2 - v0 == v2 );

 std::cout << '.';

// no multiplication scalar * vector ( only from the right side )
 assert( v1 * v2 == FT(30) );
 assert( v1 * v0 == FT(0) );
 assert( CGAL::Vector_2<R>( n1, n2) == v1 * RT(2));
 assert( CGAL::Vector_2<R>( n5, n6) == v2 * RT(3));
 assert( v2 / RT(3) == CGAL::Vector_2<R>( RT(1), -n4) );
 assert( (v2 * RT(3)) / RT(3) == v2 );
 assert( (v2 / RT(3)) * RT(3) == v2 );

#ifdef CGAL_VECTOR_MULTIPLICATION_FROM_LEFT
 assert( CGAL::Vector_2<R>( n1, n2) == RT(2) * v1 );
 assert( (CGAL::Vector_2<R>&)(RT(2) * v1) == CGAL::Vector_2<R>( n1, n2) );
 assert( CGAL::Vector_2<R>( n5, n6) == RT(3) * v2 );
#endif // CGAL_VECTOR_MULTIPLICATION_FROM_LEFT

 assert( (v4 / (FT(n1)/FT(n3))) == v10 );
 assert( (v4 * (FT(n1)/FT(n3))) == v11 );
 assert( (v4 / (FT(n3)/FT(n1))) == v11 );
 assert( (v4 * (FT(n3)/FT(n1))) == v10 );

 std::cout << '.';

 assert( v2.cartesian(0) == v2[0] );
 assert( v2.cartesian(1) == v2[1] );

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

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_VECTOR_2_C
