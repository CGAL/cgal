// Tests that all mixed operators on Quotient are defined.
// Tests that construction from double is fine.
// Sylvain Pion

#include <CGAL/Quotient.h>
#include <CGAL/use.h>
#include <limits>
#include <cassert>
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

// template class CGAL::Quotient<CGAL::Gmpz>;

template < typename T >
void test_double_ctor()
{
  typedef CGAL::Quotient<T> Qt;
  double d1 = 1e100;
  double d2 = 0.5;
  Qt q1 = d2;
  Qt q2 = d1;
  assert(q2.numerator() == T(d1));
  assert(q2.denominator() == T(1));
  assert(q1.numerator() == 1);
  assert(q1.denominator() == 2);
  q1 = d1;
  q2 = d2;
  assert(q1.numerator() == T(d1));
  assert(q1.denominator() == T(1));
  assert(q2.numerator() == 1);
  assert(q2.denominator() == 2);

  // though not exact, implicit operations +,-,* with double
  // should be the exactly the same as via explicit ctor!
  assert( (q1+0.3) == (q1+Qt(0.3)));
  assert( (q1-0.3) == (q1-Qt(0.3)));
  assert( (q1*0.3) == (q1*Qt(0.3)));
  assert( (q1/0.3) == (q1/Qt(0.3)));

  // Construction from a subnormal double used to fail with a postcondition
  // violation: no finite double can hold the denominator of denorm_min.
  if ((std::numeric_limits<double>::has_denorm)() != std::denorm_absent)
  {
    const double d = (std::numeric_limits<double>::denorm_min)();
    Qt q = d;
    assert(q.numerator() != T(0));
    T num = q.numerator();
    T den = q.denominator();
    // denorm_min == 2^-1074, hence den == num * 2^1074.
    for (int i = 0; i < 1074; ++i) num *= 2;
    assert(num == den);

    Qt qneg = -d;
    assert(qneg.numerator() == -q.numerator());
    assert(qneg.denominator() == q.denominator());
  }
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
  bool b; // avoid clang warning: equality comparison result unused [-Wunused-comparison]
  b = q<q; b = q<r; b = r<q; b = q<1; b = 1<q;
  b = q>q; b = q>r; b = r>q; b = q>1; b = 1>q;

  b = q<=q; b = q<=r; b = r<=q; b = q<=1; b = 1<=q;
  b = q>=q; b = q>=r; b = r>=q; b = q>=1; b = 1>=q;
  b = q==q; b = q==r; b = r==q; b = q==1; b = 1==q;
  b = q!=q; b = q!=r; b = r!=q; b = q!=1; b = 1!=q;
  CGAL_USE(b);
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
