// TODO : We should have a concept checker for the NT family.
// Maybe it would check too much functionalities for being useful in
// Cartesian<> and co, but it would be good for the test-suite.
// Maybe with some basic run-time testing.

// So at the moment, let me just check a few operations.

#include <CGAL/Cartesian.h>
#include <CGAL/Number_type_checker.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

typedef CGAL::Quotient<CGAL::MP_Float>                   NT0;

struct my_cmp
{
  bool operator()(const double &a, const NT0 &b) const { return NT0(a) == b; }
};

typedef CGAL::Number_type_checker<double, NT0, my_cmp>   NT;
typedef CGAL::Cartesian<NT>                              K;

int main()
{
  K::Point_2 p(0, 1);
  K::Point_2 q(2, 4);
  K::Point_2 r(6, 0);

  CGAL::orientation(p, q, r);

  return 0;
}
