// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Hervé Brönnimann

#ifndef CGAL__TEST_FCT_VECTOR_D_C
#define CGAL__TEST_FCT_VECTOR_D_C

#include <CGAL/_test_fct_vector_d.h>

template <class R>
bool
_test_fct_vector_d(const R& )
{
 std::cout << "Testing functions Vector_d" ;

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
 RT n12( 54 );


 CGAL::Vector_d<R>  v0(3, CGAL::NULL_VECTOR);  // ( 0, 0, 0)
 CGAL::Vector_d<R>  v1(3, n1, n2, n3, n4);    // ( 6,-2, 3)
 CGAL::Vector_d<R>  v2(3, n5, n6, n7, n8);    // ( 3,-6,12)
 CGAL::Vector_d<R>  v3(3, n5, n10, n9);       // ( 9,-8,15)
 CGAL::Vector_d<R>  v4(3, n8, -n2, -n5);      // ( 3, 4,-9)
 CGAL::Vector_d<R> mv4(3, -n8, n2, n5);       // (-3,-4, 9)

 CGAL::Vector_d<R>  v10(3, n8,-n2,-n5, n4);  // (1.5,2,-4.5)
 CGAL::Vector_d<R>  v11(3, -n6, n11,-n12, n8);// ( 6, 8, -18)

 assert( v1 + v2 == v3 );
 assert( v1 - v2 == v4 );
 assert( v3 - v1 == v2 );
 assert( v3 - v2 == v1 );
 assert( v4 + v2 == v1 );
 assert( v4 + v2 == v3 - v2);

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

 assert( v1 * v2 == FT(66) );
 assert( v1 * v0 == FT(0) );
 assert( CGAL::Vector_d<R>(3, n1, n2, n3) == v1 * RT(2));
 assert( CGAL::Vector_d<R>(3, n5, n6, n7) == v2 * RT(3));
 assert( v2 / RT(3) == CGAL::Vector_d<R>(3,  RT(1), -n4, -n2) );
 assert( (v2 * RT(3)) / RT(3) == v2 );
 assert( (v2 / RT(3)) * RT(3) == v2 );

 assert( (v4 / (FT(n1)/FT(n3))) == v10 );
 assert( (v4 * (FT(n1)/FT(n3))) == v11 );
 assert( (v4 / (FT(n3)/FT(n1))) == v11 );
 assert( (v4 * (FT(n3)/FT(n1))) == v10 );

 std::cout << '.';

 assert( v2.cartesian(0) == v2[0] );
 assert( v2.cartesian(1) == v2[1] );
 assert( v2.cartesian(2) == v2[2] );


 CGAL::Point_d<R> p0(3, CGAL::ORIGIN);
 CGAL::Point_d<R> p1 = CGAL::ORIGIN + v1;
 CGAL::Point_d<R> p2 = CGAL::ORIGIN + v2;
 CGAL::Point_d<R> p3 = CGAL::ORIGIN + v3;

 assert( CGAL::ORIGIN + v2 == CGAL::Point_d<R>(3,  n5, n6, n7, n8) );
 assert( CGAL::ORIGIN - v2 == CGAL::Point_d<R>(3,  -n5, -n6, -n7, n8) );
 assert( p1 - p1 == v0 );
 assert( p1 - p0 == p1 - CGAL::ORIGIN);
 assert( p1 - p2 == v4 );
 assert( p2 + v4 == p1 );
 assert( p3 - v1 == p2 );
 assert( p3 - p1 == v2 );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_VECTOR_D_C
