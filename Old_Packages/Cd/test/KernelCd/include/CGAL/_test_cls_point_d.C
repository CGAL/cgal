#ifndef CGAL__TEST_CLS_POINT_D_H
#define CGAL__TEST_CLS_POINT_D_H

#include <CGAL/_test_cls_point_d.h>

template <class R>
bool
_test_cls_point_d( const R& )
{
  std::cout << "Testing class Point_d" ;

  double coord1[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 1.0};
  double coord2[6] = {0.0, 2.0, 4.0, 6.0, 8.0, 2.0};
  CGAL::Point_d<R> p (5, coord1, coord1+6); 
  CGAL::Point_d<R> pp (5, coord1, coord1+5);
  CGAL::Point_d<R> q (5, coord2, coord2+6);
  CGAL::Point_d<R> s (p);   // constructors
  CGAL::Point_d<R> t;
  t = s;                    // assignment

  std::cout << '.';
  
  assert ( p == pp);
  assert ( p == q);
  assert ( p == s);         // equality test
  assert ( p == t);

  std::cout << '.';
  
  for (int i=0; i<5; ++i) 
  {
    assert (q.cartesian(i) == (double)i);   // cartesian method
    assert (q[i] == (double)i);             // operator[]
    assert (p.homogeneous(i) == (double)i);
  }
  assert (s.dimension() == 5);              // dimension

  std::cout << "done" << std::endl;
  return true;
}  
 
#endif // CGAL__TEST_CLS_POINT_D_H
