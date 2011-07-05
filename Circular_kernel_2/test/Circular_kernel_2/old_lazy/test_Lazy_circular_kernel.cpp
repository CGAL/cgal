#include <CGAL/Cartesian.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_circular_kernel_2.h>
#include <CGAL/intersections.h>
#include <iostream>

typedef CGAL::Quotient<CGAL::MP_Float>                       NT3;
typedef CGAL::Cartesian<NT3>                                 Linear_k3;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT3>          Algebraic_k3;
typedef CGAL::Circular_kernel_2<Linear_k3, Algebraic_k3>     CK1_2;
  
typedef CGAL::Interval_nt_advanced                           NT4;
typedef CGAL::Cartesian<NT4>                                 Linear_k4;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT4>          Algebraic_k4;
typedef CGAL::Circular_kernel_2<Linear_k4,Algebraic_k4>      CK2_2;
typedef CGAL::Lazy_circular_kernel_2<CK1_2,CK2_2>            CK2;

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
