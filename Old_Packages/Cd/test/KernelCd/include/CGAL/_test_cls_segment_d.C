// revision      : 2.2.2
// revision_date : 28 Sep 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>

#ifndef CGAL__TEST_CLS_SEGMENT_D_C
#define CGAL__TEST_CLS_SEGMENT_D_C

#include <CGAL/_test_cls_segment_d.h>

template <class R>
bool
_test_cls_segment_d(const R& )
{
 std::cout << "Testing class Segment_d" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Segment_d is;
 CGAL::Segment_d<R>  s1(is);

 RT  n1 =  7;
 RT  n2 = 21;
 RT  n3 = 14;
 RT  n4 =-10;
 RT  n5 =  5;
 RT  n6 = 20;
 RT  n7 = -1;
 RT  n8 =  3;

 CGAL::Point_d<R> p1(3, n1, n2, n3, n1);
 CGAL::Point_d<R> p2(3, n4, n5, n6, n5);
 CGAL::Point_d<R> p3(3, n2, n8, n2, n8);

 CGAL::Segment_d<R> s2( p1, p2 );
 CGAL::Segment_d<R> s3( p2, p1 );
 CGAL::Segment_d<R> s4( s2 );
 s1 = s4;

 assert( s1 == s1 );
 assert( s4 == s2 );
 assert( s1 == s4 );
 assert( s1 == s2 );
 assert( s2 != s3 );

 CGAL::Segment_d<R> s5 (p3, p3 + (p1 - p3) + (p1 - p3) );
 assert( s5.has_on( p1 ) );
 assert( s5.has_on( p3 ) );
 assert( s2.has_on( p2 ) );

 std::cout <<'.';

 assert( s5.source() == p3 );
 assert( s5.target() == p1 + (p1 - p3) );
 assert( s2.source() != s3.source() );
 assert( s2.target() != s3.target() );

 std::cout <<'.';

 assert( s2.min() == p2 );
 assert( s2.max() == p1 );
 assert( s2.min() == s3.min() );
 assert( s2.max() == s3.max() );
 assert( s5.max() != s5.min() );
 assert( s5.max() == s5.opposite().max() );
 assert( s5.vertex(0) == s5.source() );
 assert( s2.vertex(1) == s2.target() );
 assert( s2.vertex(1) == s2.min() );
 assert( s2[1] == s1[1] );
 assert( s2[1] == s3[0] );

 std::cout << '.';

 assert( s2.squared_length() == FT( RT(17) ) );
 assert( s2.direction() == CGAL::Direction_d<R>(s2.target() - s2.source() ));
 assert( s2.direction() == s3.opposite().direction() );

 assert( s1.supporting_line() == s2.supporting_line() );
 CGAL::Line_d<R> lin(p1,p2);
 assert( s2.supporting_line() == lin );

 CGAL::Segment_d<R> sdeg(p3,p3);
 assert( sdeg.is_degenerate() );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_SEGMENT_D_C
