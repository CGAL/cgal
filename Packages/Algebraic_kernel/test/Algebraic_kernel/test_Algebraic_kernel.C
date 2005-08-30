#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_2_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Quotient.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Gmpq.h>
#include <CGAL/_test_predicates.h>
#include <CGAL/_test_constructor.h>

int main()
{
  typedef CGAL::Gmpq                                           NT1;
  typedef CGAL::Algebraic_kernel_2_2<NT1>                      Algebraic_k1;
  Algebraic_k1 ak1;
  _test_solve(ak1);
  _test_sign_at(ak1);
  _test_critical_points(ak1);
  _test_compare_Root_for_circles(ak1);
  _test_constuctor(ak1);
  
  typedef CGAL::Quotient<CGAL::MP_Float>                       NT2;
  typedef CGAL::Algebraic_kernel_2_2<NT2>                      Algebraic_k2;
  Algebraic_k2 ak2;
  _test_solve(ak2);
  _test_sign_at(ak2);
  _test_critical_points(ak2);
  _test_compare_Root_for_circles(ak2);
  _test_constuctor(ak2);
  
  return 0;
}
