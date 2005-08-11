#include <CGAL/Cartesian.h>
#include <CGAL/Circular_kernel.h>
#include <CGAL/Algebraic_kernel_2_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Quotient.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Gmpq.h>
#include <CGAL/_test_circles_predicates.h>
#include <CGAL/_test_circles_constructions.h>
#include <CGAL/Filtered_hexagon_curved_kernel.h>
#include <iostream>



int main()
{
  typedef CGAL::Gmpq                                           NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_2_2<NT1>                      Algebraic_k1;
  typedef CGAL::Curved_kernel<Linear_k1,Algebraic_k1>          CK1_;
  typedef CGAL::Filtered_hexagon_curved_kernel<CK1_>           CK1;

  CK1 ck1;
  _test_circle_predicat(ck1);
  
  _test_circle_construct(ck1);

  typedef CGAL::Quotient<CGAL::MP_Float>                       NT2;
  typedef CGAL::Cartesian<NT2>                                 Linear_k2;
  typedef CGAL::Algebraic_kernel_2_2<NT2>                      Algebraic_k2;
  typedef CGAL::Curved_kernel<Linear_k2, Algebraic_k2>         CK2_;
  typedef CGAL::Filtered_hexagon_curved_kernel<CK2_>           CK2;

  
  CK2 ck2;

  CK2::Circle_2 c(CK2::Point_2(0,0),1);
  CK2::Circular_arc_2 a(c);

   std::cout<<" is y monotone"<<a.is_y_monotone()<<std::endl;

  _test_circle_predicat(ck2);
  _test_circle_construct(ck2);

  return 0;
}
