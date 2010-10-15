#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/_test_predicates.h>
#include <CGAL/_test_constructor.h>

int main()
{
  typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1>          Algebraic_k1;
  Algebraic_k1 ak1;
  _test_solve(ak1);
  _test_sign_at(ak1);
  _test_critical_points(ak1);
  _test_compare_Root_for_circles(ak1);
  _test_constuctor(ak1);

  return 0;
}
