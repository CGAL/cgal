// Tests that all mixed operators on Quotient are defined.
// Tests that construction from double is fine.
// Sylvain Pion

#include <CGAL/Quotient.h>
#include <CGAL/Testsuite/assert.h>
#include <CGAL/MP_Float.h>
#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#endif
#ifdef CGAL_USE_GMPXX
#  include <CGAL/mpz_class.h>
#endif
#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_integer.h>
#endif

template < typename T >
void test_double_ctor()
{
  typedef CGAL::Quotient<T> Qt;
  double d1 = 1e100;
  double d2 = 0.5;
  Qt q1 = d1;
  Qt q2 = d2;
  CGAL_test_assert(q1.numerator() == T(d1));
  CGAL_test_assert(q2.numerator() == 1);
  CGAL_test_assert(q2.denominator() == 2);
}

void test_comparison_operators()
{
  typedef CGAL::MP_Float     RT;
  typedef CGAL::Quotient<RT> QT;

  RT r(1);
  QT q(1);

  q+q; q+r; r+q; q+1; 1+q;
  q-q; q-r; r-q; q-1; 1-q;
  q*q; q*r; r*q; q*1; 1*q;
  q/q; q/r; r/q; q/1; 1/q;
  -q;
  q<q; q<r; r<q; q<1; 1<q;
  q>q; q>r; r>q; q>1; 1>q;
  q<=q; q<=r; r<=q; q<=1; 1<=q;
  q>=q; q>=r; r>=q; q>=1; 1>=q;
  q==q; q==r; r==q; q==1; 1==q;
  q!=q; q!=r; r!=q; q!=1; 1!=q;
}

int main()
{
  test_comparison_operators();

  // Quotient<MP_Float and Gmpzf>'s double ctor do not split the double
  // in integral parts.  So, not tested here.
#ifdef CGAL_USE_GMP
  test_double_ctor<CGAL::Gmpz>();
#endif
#ifdef CGAL_USE_GMPXX
  test_double_ctor<mpz_class>();
#endif
#ifdef CGAL_USE_LEDA
  test_double_ctor<leda_integer>();
#endif

  return 0;
}
