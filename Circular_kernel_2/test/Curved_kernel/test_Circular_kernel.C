#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_2_2.h>
#include <CGAL/Circular_kernel.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
//#include <CGAL/Gmpq.h>
#include <CGAL/_test_circles_predicates.h>
#include <CGAL/_test_circles_constructions.h>
#include <CGAL/_test_circles_extention.h>

int main()
{
  typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
  //  typedef CGAL::Gmpq                                           NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1>          Algebraic_k1;
  typedef CGAL::Circular_kernel_2<Linear_k1,Algebraic_k1>      CK1;
  CK1 ck1;
  std::cout << "Testing predicates..." << std::endl;
  _test_circle_predicat(ck1);
  std::cout << "Testing constructions..." << std::endl;
  _test_circle_construct(ck1);
  std::cout << "Testing bboxes..." << std::endl;
  _test_circle_bbox(ck1);
  std::cout << "Testing circular_arc_bboxes..." << std::endl;
  _test_circular_arc_bbox(ck1);
  std::cout << "Testing circular_arc_point_bboxes..." << std::endl;
  _test_circular_arc_point_bbox(ck1);
  std::cout << "Testing has_on..." << std::endl;
  _test_has_on(ck1);
/*
  typedef CGAL::Quotient<CGAL::MP_Float>                       NT2;
  typedef CGAL::Cartesian<NT2>                                 Linear_k2;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT2>          Algebraic_k2;
  typedef CGAL::Circular_kernel_2<Linear_k2, Algebraic_k2>     CK2;
  
  CK2 ck2;
  _test_circle_predicat(ck2);
  _test_circle_construct(ck2);
  _test_circle_bbox(ck2);
  _test_circular_arc_bbox(ck2);
  _test_has_on(ck2);*/
  return 0;
}
