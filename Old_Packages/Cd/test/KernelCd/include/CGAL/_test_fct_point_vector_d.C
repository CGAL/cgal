#ifndef CGAL__TEST_FCT_POINT_VECTOR_D_C
#define CGAL__TEST_FCT_POINT_VECTOR_D_C

#include <CGAL/_test_fct_point_vector_d.h>

template <class R>
bool
_test_fct_point_vector_d(const R& )
{
 std::cout << "Testing functions Point_d Vector_d" ;

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

 CGAL::Vector_d<R>  v0(3, CGAL::NULL_VECTOR);
 CGAL::Vector_d<R>  v1(3, n1, n2, n3, n4);
 CGAL::Vector_d<R>  v2(3, n5, n6, n7, n8);
 CGAL::Vector_d<R>  v3(3, n5, n10, n9);
 CGAL::Vector_d<R>  v4(3, n8, -n2, -n5);

 std::cout << '.';

 CGAL::Point_d<R> p0(3, CGAL::ORIGIN);
 CGAL::Point_d<R> p1 = CGAL::ORIGIN + v1;
 CGAL::Point_d<R> p2 = CGAL::ORIGIN + v2;
 CGAL::Point_d<R> p3 = CGAL::ORIGIN + v3;

 assert( CGAL::ORIGIN + v2 == CGAL::Point_d<R>(3, n5, n6, n7, n8) );
 assert( CGAL::ORIGIN - v2 == CGAL::Point_d<R>(3, -n5, -n6, -n7, n8) );
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
#endif // CGAL__TEST_FCT_POINT_VECTOR_D_C
