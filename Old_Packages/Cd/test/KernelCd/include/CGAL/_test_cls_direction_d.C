#ifndef CGAL__TEST_CLS_DIRECTION_D_C
#define CGAL__TEST_CLS_DIRECTION_D_C

#include <CGAL/_test_cls_direction_d.h>

template <class R>
bool
_test_cls_direction_d(const R& )
{
 std::cout << "Testing class Direction_d" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Direction_d  id;

 CGAL::Direction_d<R> d0;
 CGAL::Direction_d<R> d1(id);

 std::cout << '.';
 RT   n0 = 10;
 RT  n1 = 8;
 RT  n2 = 4;
 RT  n3 = 2;

 CGAL::Vector_d<R>  v(3, n1, n2, n3);
 CGAL::Direction_d<R> d2(v);
 CGAL::Direction_d<R> d3(3, n0, n1, n2);
 CGAL::Direction_d<R> d4( d3 );
 CGAL::Direction_d<R> d5 = d3;

 assert( d3 == d3 );
 assert( d3 == d4 );
 assert( d5 == d3 );
 assert( d2 != d3 );
 assert( d3 != d2 );

 std::cout << '.';
 CGAL::Vector_d<R> vv = d2.to_vector();
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

#endif // CGAL__TEST_CLS_DIRECTION_D_C
