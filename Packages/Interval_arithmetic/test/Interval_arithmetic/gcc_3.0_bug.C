/*
  This test file has been distilled from a bug-report from Matthias.
  Probably a GCC 3.0.4 code generation bug.
  It shows up only with -O -DNDEBUG.
  It goes away when :
  - operator-(Interval) is not inline
  - copy ctor of Interval_base is written explicitly.

> Have you made some progress for this bug ?
> If you could try to keep the bug without having the test program depend on
> LEDA, I could have a look at it, but currently I can't even compile it...
> -- 
> Sylvain

I've changed one of the examples from the Triangulation_3 examples a bit
and it shows the same behavior.

When I compile with -DNDEBUG it loops; removing -DNDEBUG from the
makefile makes it work. I attach the makefile and the example.

Matthias

*/

#ifndef NDEBUG
#define NDEBUG
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Interval_arithmetic.h>

#include <iostream>

typedef CGAL::Interval_nt_advanced NT;
typedef CGAL::Simple_cartesian<NT> K;

int main()
{
  std::cout.precision(20);
  CGAL::Bounded_side bs;
  K::Point_3 p0(0, 0, 0);
  K::Point_3 p1(1, 0, 0);
  K::Point_3 p2(0, 1, 0);

 {
  CGAL::Protect_FPU_rounding<false> p;
  bs = CGAL::coplanar_side_of_bounded_circle(p0, p1, p2, p0);
 }

  if (bs != CGAL::ON_BOUNDARY) {
    std::cout << "BUG !  " << (int) bs << std::endl;
    abort();
  }

  std::cout << p0 << std::endl;
  std::cout << p1 << std::endl;
  std::cout << p2 << std::endl;

  return 0;
}
