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
// file          : _test_fct_direction_2.C
// revision      : 2.1
// revision_date : 05 Aug 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL__TEST_FCT_DIRECTION_2_C
#define CGAL__TEST_FCT_DIRECTION_2_C
#ifndef CGAL__TEST_FCT_DIRECTION_2_H
#include <CGAL/_test_fct_direction_2.h>
#endif // CGAL__TEST_FCT_DIRECTION_2_H

template <class R>
bool
_test_fct_direction_2(const R& )
{
 std::cout << "Testing functions Direction_2" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 RT n0 =  0;
 RT n1 =  1;
 RT n2 =  2;
 RT n3 =  3;
 RT n4 =  4;
 RT n5 =  5;
 RT n6 =  6;
 RT n8 =  8;

 CGAL::Direction_2<R> d[8];
 d[0] = CGAL::Direction_2<R>( n5, n0 );
 d[1] = CGAL::Direction_2<R>( n8, n3 );
 d[2] = CGAL::Direction_2<R>( n0, n4 );
 d[3] = CGAL::Direction_2<R>(-n4, n3 );
 d[4] = CGAL::Direction_2<R>(-n2, n0 );
 d[5] = CGAL::Direction_2<R>(-n4,-n6 );
 d[6] = CGAL::Direction_2<R>( n0,-n1 );
 d[7] = CGAL::Direction_2<R>( n4,-n5 );

 std::cout << '.';

 assert( d[0] >= d[0] );
 assert( d[0] <= d[0] );

 std::cout << '.';

 int i;
 int j;

 for ( i = 1; i <= 7; i++ )
 {
    // std::cout << std::endl;
    for ( j = 0; j+i <= 7; j++)
    {
        // std::cout << '('  << j << ',' << j+i << ')' ;
        assert( d[j] <= d[j+i] );
        assert( d[j] <  d[j+i] );
    }
    // std::cout << ' ';
    for (      ; j <= 7; j++)
    {
        // std::cout << '('  << j << ',' << (j+i)%8 << ')' ;
        assert( d[j] >= d[(j+i)%8] );
        assert( d[j] >  d[(j+i)%8] );
    }
    assert( d[i] >= d[i] );
    assert( d[i] <= d[i] );
 }

 std::cout << '.';

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_DIRECTION_2_C
