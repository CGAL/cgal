#include <utility>
namespace CGAL {
class MP_Float;
template < typename T > class Root_of_2;
template < typename T > class Lazy_exact_nt;
template < typename T >
std::pair<double,double> to_interval(const Root_of_2<T>&);

#if 0
template < typename T >
void operator-(Lazy_exact_nt<Root_of_2<T> > a, Lazy_exact_nt< Root_of_2<T> > b) { f(a); }
void operator-(Lazy_exact_nt<Root_of_2<MP_Float> > a, Lazy_exact_nt< Root_of_2<MP_Float> > b);

template < typename T1, typename T2 >
struct Binary_operator_result;

template < typename T1, typename T2 >
struct Binary_operator_result <Root_of_2<T1>, Root_of_2<T2> >;
#endif
}

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Circular_kernel.h>
#include <CGAL/Algebraic_kernel_2_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Quotient.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Gmpq.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Lazy_exact_nt.h>
#include <CGAL/Lazy_curved_kernel.h>
#include <CGAL/intersections.h>
#include <iostream>



//int main()
//{
  typedef CGAL::Gmpq                                           NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_2_2<NT1>                      Algebraic_k1;
  typedef CGAL::Curved_kernel<Linear_k1,Algebraic_k1>          CK1_1;
  
  typedef CGAL::Interval_nt<>                                  NT2;
  typedef CGAL::Cartesian<NT2>                                 Linear_k2;
  typedef CGAL::Algebraic_kernel_2_2<NT2>                      Algebraic_k2;
  typedef CGAL::Curved_kernel<Linear_k2,Algebraic_k2>          CK2_1;
  typedef CGAL::Lazy_curved_kernel<CK1_1,CK2_1>                  CK1;
  
// CK1 ck1;
//  _test_circle_predicat(ck1);
//  _test_circle_construct(ck1);

  typedef CGAL::Quotient<CGAL::MP_Float>                       NT3;
  typedef CGAL::Cartesian<NT3>                                 Linear_k3;
  typedef CGAL::Algebraic_kernel_2_2<NT3>                      Algebraic_k3;
  typedef CGAL::Curved_kernel<Linear_k3, Algebraic_k3>         CK1_2;
  
  //typedef CGAL::Interval_nt<>                                  NT4;
  typedef CGAL::Interval_nt_advanced                           NT4;
  typedef CGAL::Cartesian<NT4>                                 Linear_k4;
  typedef CGAL::Algebraic_kernel_2_2<NT4>                      Algebraic_k4;
  typedef CGAL::Curved_kernel<Linear_k4,Algebraic_k4>          CK2_2;
  typedef CGAL::Lazy_curved_kernel<CK1_2,CK2_2>                  CK2;

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
