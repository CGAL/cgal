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
// file          : _test_cls_vector_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_CLS_VECTOR_2_C
#define CGAL__TEST_CLS_VECTOR_2_C

#ifndef CGAL__TEST_CLS_VECTOR_2_H
#include <CGAL/_test_cls_vector_2.h>
#endif // CGAL__TEST_CLS_VECTOR_2_H

template <class R>
bool
_test_cls_vector_2(const R& )
{
 std::cout << "Testing class Vector_2" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Vector_2       iv;
 CGAL::Vector_2<R>  v1;
 CGAL::Vector_2<R>  v2(iv);
 CGAL::Vector_2<R>  v0(CGAL::NULL_VECTOR);

 RT  n1( 12 );
 RT  n2( -4 );
 RT  n3(  6 );
 RT  n4(  2 );

 CGAL::Vector_2<R>  v3(n1, n2 );       // ( 12, -4)
 CGAL::Vector_2<R>  v4(n1, n2, n4);    // (  6, -2)
 CGAL::Vector_2<R>  v5(n1, n2, n4);    // (  6, -2)
 CGAL::Vector_2<R>  v6( v5 );
                   v1 = v4;
 CGAL::Vector_2<R>  v7(-n1, -n2, -n4); // (  6, -2)

 std::cout << '.';

 assert( v5 == v7 );
 assert( v4 == v5 );
 assert( v5 == v6 );
 assert( v4 == v6 );
 assert( v1 == v6 );
 assert( v0 == CGAL::NULL_VECTOR);

 assert( v3 != v4 );
 assert( v0 != v1 );
 assert( v1 != CGAL::NULL_VECTOR);

 assert( v3.hx() == n1 );   // don't replace v3
 assert( v3.hy() == n2 );

 assert( FT( v5.hx()) / FT(v5.hw()) == FT( n1) / FT( n4) );
 assert( FT( v5.hy()) / FT(v5.hw()) == FT( n2) / FT( n4) );

 assert( v5.x() == FT( n1) / FT( n4) );
 assert( v5.y() == FT( n2) / FT( n4) );

 std::cout << '.';

 assert( v3.homogeneous(0) == v3.hx() );  // don't replace v3
 assert( v3.homogeneous(1) == v3.hy() );
 assert( v3.homogeneous(2) == v3.hw() );
 assert( v6.cartesian(0) == v6.x() );
 assert( v6.cartesian(1) == v6.y() );

 std::cout << '.';

 assert( v0.dimension() == 2 );
 assert( v4.homogeneous( v4.dimension() ) == v4.hw() );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_VECTOR_2_C
