
#include <CGAL/CORE_Expr.h>

typedef CORE::Expr Real;

int main()
{
  Real r(3.14);

  CGAL::is_zero(r);

  CGAL::is_one(r);

  r = CGAL::sqrt(r);


  CGAL::is_zero(r);

  CGAL::is_one(r);

  r = r * r;

  CGAL::is_zero(r);

  CGAL::is_one(r);

  r = r - r;

  CGAL::is_zero(r);

  return 0;
}
