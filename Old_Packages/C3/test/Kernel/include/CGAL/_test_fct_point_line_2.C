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
// file          : _test_fct_point_line_2.C
// revision      : 2.0.5
// revision_date : 24 Mar 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL__TEST_FCT_POINT_2_LINE_2_C
#define CGAL__TEST_FCT_POINT_2_LINE_2_C
#ifndef CGAL__TEST_FCT_POINT_LINE_2_H
#include <CGAL/_test_fct_point_line_2.h>
#endif // CGAL__TEST_FCT_POINT_LINE_2_H

template <class R>
bool
_test_fct_point_line_2(const R& )
{
 cout << "Testing functions Point_2 Line_2" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 RT n0 =  0;
 RT n1 =  1;
 RT n2 =  2;
 RT n3 =  3;
 RT n4 =  4;
 RT n5 =  5;
 RT n6 =  6;
 RT n7 =  7;
 RT n8 =  8;
 RT n9 =  9;
 RT n10= 10;
 RT n12= 12;
 RT n20= 20;

 CGAL::Point_2<R>  p1( n2, n8, n2);   // ( 1, 4)
 CGAL::Point_2<R>  p2( n4, n4, n2);   // ( 2, 2)
 CGAL::Point_2<R>  p3(n12,n10, n2);   // ( 6, 5)
 CGAL::Point_2<R>  p4( n7, n3);       // ( 7, 3)
 CGAL::Point_2<R>  p5( n9, n0, n3);   // ( 3, 0)
 CGAL::Point_2<R>  p6(n12,n20, n4);   // ( 3, 5)
 CGAL::Point_2<R>  p7(n12, n3, n3);   // ( 4, 1)
 CGAL::Point_2<R>  p8( n8, n4);       // ( 8, 4)
 CGAL::Point_2<R>  p9(-n4, n8, n4);   // (-1, 2)
 CGAL::Point_2<R>  p10(n10, n4, n2);  // ( 5, 2)

 CGAL::Point_2<R>  p11( n1, -n5,-n1); // (-1, 5)
 CGAL::Point_2<R>  p0( CGAL::ORIGIN ); // ( 0, 0)
 CGAL::Point_2<R>  p13( n8, n8, n4);  // ( 4, 4)
 CGAL::Point_2<R>  p14( n5,-n1);      // ( 5,-1)
 CGAL::Point_2<R>  p16(n12, n9, n3);  // ( 4, 3)
 CGAL::Point_2<R>  p17( n0,n1);       // ( 0, 1)

 CGAL::Line_2<R>   l14( p1, p4 );
 CGAL::Line_2<R>   l23( p2, p3 );
 CGAL::Line_2<R>   l67( p6, p7 );
 CGAL::Line_2<R>   l58( p5, p8 );
 CGAL::Line_2<R>   l1617( p16, p17);
 CGAL::Line_2<R>   l1114( p11, p14);
 CGAL::Line_2<R>   l1716( p17, p16);
 CGAL::Line_2<R>   l1411( p14, p11);
 CGAL::Line_2<R>   l013( p0, p13 );
 CGAL::Line_2<R>   l910( p9, p10 );

 cout << '.';

 assert( CGAL::compare_x( p16, p14 ) == CGAL::SMALLER );

 assert( CGAL::compare_x( p9, l14, l23) == CGAL::SMALLER );
 assert( CGAL::compare_x( p8, l14, l23) == CGAL::LARGER );
 assert( CGAL::compare_x( p2, l1617, l910) == CGAL::EQUAL );
 assert( CGAL::compare_x( p2, l1716, l910) == CGAL::EQUAL );
 assert( CGAL::compare_x( p2, l1114, l013) == CGAL::EQUAL );

 assert( CGAL::compare_y( p6, l14, l23 ) == CGAL::LARGER );
 assert( CGAL::compare_y( p9, l14, l23 ) == CGAL::SMALLER );
 assert( CGAL::compare_x( p2, l1411, l013) == CGAL::EQUAL );
 assert( CGAL::compare_x( p2, l1716, l013) == CGAL::EQUAL );

 cout << '.';

 assert( CGAL::compare_x( l14, l23, l58 ) == CGAL::SMALLER);
 assert( CGAL::compare_x( l14, l58, l23 ) == CGAL::LARGER);
 assert( CGAL::compare_x( l14, l58, l58 ) == CGAL::EQUAL);
 assert( CGAL::compare_x( l1114, l013, l910 ) == CGAL::EQUAL);
 assert( CGAL::compare_x( l1617, l910, l013 ) == CGAL::EQUAL);

 assert( CGAL::compare_y( l14, l58, l23 ) == CGAL::SMALLER);
 assert( CGAL::compare_y( l14, l23, l58 ) == CGAL::LARGER);
 assert( CGAL::compare_y( l14, l58, l58 ) == CGAL::EQUAL);
 assert( CGAL::compare_y( l1114, l013, l910 ) == CGAL::EQUAL);
 assert( CGAL::compare_y( l1617, l910, l013 ) == CGAL::EQUAL);

 assert( CGAL::compare_x( l14, l23, l67, l58 ) == CGAL::SMALLER);
 assert( CGAL::compare_x( l67, l58, l23, l14 ) == CGAL::LARGER);
 assert( CGAL::compare_x( l1114, l1617, l910, l013 ) == CGAL::EQUAL);

 assert( CGAL::compare_y( l14, l23, l67, l58 ) == CGAL::LARGER);
 assert( CGAL::compare_y( l67, l58, l23, l14 ) == CGAL::SMALLER);
 assert( CGAL::compare_y( l1114, l1617, l910, l013 ) == CGAL::EQUAL);

 cout << '.';

 assert( CGAL::compare_y_at_x( p6, l23 ) == CGAL::LARGER );
 assert( CGAL::compare_y_at_x( p6, l23.opposite() ) == CGAL::LARGER );
 assert( CGAL::compare_y_at_x( p10, l23 ) == CGAL::SMALLER );
 assert( CGAL::compare_y_at_x( p9, l23 ) == CGAL::LARGER );
 assert( CGAL::compare_y_at_x( p17, l910 ) == CGAL::SMALLER );
 assert( CGAL::compare_y_at_x( p0, l23 ) == CGAL::SMALLER );
 assert( CGAL::compare_y_at_x( p8, l58 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( p2, l1617 ) == CGAL::EQUAL );

 assert( CGAL::compare_y_at_x( l14, l23, l58 ) == CGAL::LARGER );
 assert( CGAL::compare_y_at_x( l67, l58, l23 ) == CGAL::SMALLER );
 assert( CGAL::compare_y_at_x( l910, l1716, l1114) == CGAL::EQUAL);

 assert( CGAL::compare_y_at_x( l14, l23, l58, l67 ) == CGAL::SMALLER );
 assert( CGAL::compare_y_at_x( l14, l23, l67, l58 ) == CGAL::LARGER );
 assert( CGAL::compare_y_at_x( l14, l23, l1411, l1114 ) == CGAL::EQUAL );
 assert( CGAL::compare_y_at_x( l1617, l013, l910, l67 ) == CGAL::SMALLER);
 assert( CGAL::compare_y_at_x( l1617, l013, l67, l910 ) == CGAL::LARGER);
 assert( CGAL::compare_y_at_x( l1617, l013, l1114, l910 ) == CGAL::EQUAL);
 assert( CGAL::compare_y_at_x( l1617, l013, l910, l1114 ) == CGAL::EQUAL);

 cout << "done" << endl;
 return true;
}
#endif // CGAL__TEST_FCT_POINT_2_LINE_2_C
