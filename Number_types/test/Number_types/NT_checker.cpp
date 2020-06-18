// TODO : We should have a concept checker for the NT family.
// Maybe it would check too much functionalities for being useful in
// Cartesian<> and co, but it would be good for the test-suite.
// Maybe with some basic run-time testing.

// So at the moment, let me just check a few operations.

#include <CGAL/Cartesian.h>
#include <CGAL/Number_type_checker.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>

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

#ifdef CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
  { // Integer
      typedef CGAL::Arithmetic_kernel::Integer Integer;
      typedef CGAL::Number_type_checker<Integer,Integer> NT;
      typedef CGAL::Euclidean_ring_tag Tag;
      typedef CGAL::Tag_true Is_exact;
      CGAL::test_algebraic_structure<NT,Tag, Is_exact>();
  }
  {// Rational
      typedef CGAL::Arithmetic_kernel::Rational Rational;
      typedef CGAL::Number_type_checker<Rational,Rational> NT;
      typedef CGAL::Field_tag Tag;
      typedef CGAL::Tag_true Is_exact;
      CGAL::test_algebraic_structure<NT,Tag, Is_exact>();
  }
#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
  return 0;
}
