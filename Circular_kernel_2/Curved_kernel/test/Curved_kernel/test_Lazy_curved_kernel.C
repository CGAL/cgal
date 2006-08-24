#include <CGAL/Cartesian.h>
#include <CGAL/Circular_kernel.h>
#include <CGAL/Algebraic_kernel_2_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
//#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_curved_kernel.h>
#include <CGAL/intersections.h>
#include <iostream>



//int main()
//{
  typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
//  typedef CGAL::Gmpq                                           NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1>          Algebraic_k1;
  typedef CGAL::Circular_kernel_2<Linear_k1,Algebraic_k1>      CK1_1;
  
  typedef CGAL::Interval_nt<>                                  NT2;
  typedef CGAL::Cartesian<NT2>                                 Linear_k2;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT2>          Algebraic_k2;
  typedef CGAL::Circular_kernel_2<Linear_k2,Algebraic_k2>      CK2_1;
  typedef CGAL::Lazy_curved_kernel<CK1_1,CK2_1>                CK1;
  
// CK1 ck1;
//  _test_circle_predicat(ck1);
//  _test_circle_construct(ck1);

  typedef CGAL::Quotient<CGAL::MP_Float>                       NT3;
  typedef CGAL::Cartesian<NT3>                                 Linear_k3;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT3>          Algebraic_k3;
  typedef CGAL::Circular_kernel_2<Linear_k3, Algebraic_k3>     CK1_2;
  
  //typedef CGAL::Interval_nt<>                                  NT4;
  typedef CGAL::Interval_nt_advanced                           NT4;
  typedef CGAL::Cartesian<NT4>                                 Linear_k4;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT4>          Algebraic_k4;
  typedef CGAL::Circular_kernel_2<Linear_k4,Algebraic_k4>      CK2_2;
  typedef CGAL::Lazy_curved_kernel<CK1_2,CK2_2>                CK2;

  CK2 ck2;

#include <CGAL/_test_circles_predicates.h>
#include <CGAL/_test_circles_constructions.h>
#include <CGAL/_test_circles_extention.h>
  
int main() {

  _test_circle_predicat(ck2);
  _test_circle_construct(ck2);
  _test_circle_bbox(ck2);
  _test_circular_arc_bbox(ck2);
  _test_has_on(ck2);

  return 0;
}
