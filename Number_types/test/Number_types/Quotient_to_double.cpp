// Test for the generic Real_embeddable_traits<Quotient<NT> >::To_double,
// covering issue #1053 and issue #1815.
//
// #1053: the old implementation divided in NT, but NT is only required to be a
//        model of IntegralDomainWithoutDivision, so it need not have operator/.
// #1815: the old implementation tested is_finite() on the NT values rather than
//        on the converted doubles.  For an integer type is_finite() is always
//        true, so the overflow handling never fired and the result was NaN.
//
// The conversion is to_double(num) / to_double(den), and out of range values are
// reported the way IEEE reports them.  Note that this conversion is not
// monotonic: the numerator and the denominator are rounded independently, so two
// close rationals can be converted in the opposite order.  That is issue #8966,
// which is about Quotient<Gmpzf> and its own specialization of To_double, and it
// is not addressed here.

#include <cassert>
#include <cmath>
#include <iostream>

#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif

template <typename NT>
void test_basic()
{
  typedef CGAL::Quotient<NT> Q;

  assert(CGAL::to_double(Q(NT(0), NT(1))) == 0.0);
  assert(CGAL::to_double(Q(NT(42), NT(1))) == 42.0);
  assert(CGAL::to_double(Q(NT(1), NT(2))) == 0.5);
  assert(CGAL::to_double(Q(NT(-1), NT(2))) == -0.5);
  assert(CGAL::to_double(Q(NT(3), NT(4))) == 0.75);

  const double d = CGAL::to_double(Q(NT(-7), NT(3)));
  assert(d > -2.4 && d < -2.3);

  std::cout << "  basic conversions: OK" << std::endl;
}

#ifdef CGAL_USE_GMP
// Build 10^e as a Gmpz using multiplication only.
CGAL::Gmpz pow10(int e)
{
  CGAL::Gmpz r(1);
  for (int i = 0; i < e; ++i)
    r = r * CGAL::Gmpz(10);
  return r;
}

// Quotient<Gmpz> uses the generic To_double above, so it exercises the code
// this test is about.
void test_out_of_range()
{
  typedef CGAL::Quotient<CGAL::Gmpz> Q;

  const CGAL::Gmpz huge = pow10(400);          // far beyond the double range
  assert(!CGAL::is_finite(CGAL::to_double(huge)));

  // Only the denominator is out of range: the quotient underflows to zero.
  {
    const double d = CGAL::to_double(Q(CGAL::Gmpz(1), huge));
    assert(d == 0.0);
  }

  // Only the numerator is out of range: the quotient overflows to infinity.
  {
    const double d = CGAL::to_double(Q(huge, CGAL::Gmpz(2)));
    assert(!CGAL::is_finite(d));
    assert(d > 0.0);
  }

  // Negative numerator, same thing with the opposite sign.
  {
    const double d = CGAL::to_double(Q(-huge, CGAL::Gmpz(2)));
    assert(!CGAL::is_finite(d));
    assert(d < 0.0);
  }

  // Both out of range.  Nothing can be recovered from the two converted
  // values alone, and the result is NaN rather than a made up finite value.
  // Before the fix this case was reached through a division in NT.
  {
    const double d = CGAL::to_double(Q(huge, huge + CGAL::Gmpz(1)));
    assert(std::isnan(d));
  }

  std::cout << "  out of range conversions: OK" << std::endl;
}

// Values whose numerator and denominator both stay inside the double range are
// converted through the plain division, with no NT division involved.
void test_large_but_representable()
{
  typedef CGAL::Quotient<CGAL::Gmpz> Q;

  const CGAL::Gmpz a = pow10(40);
  const CGAL::Gmpz b = pow10(38);

  const double d = CGAL::to_double(Q(a, b));   // 10^2
  assert(CGAL::is_finite(d));
  assert(d > 99.0 && d < 101.0);

  const double e = CGAL::to_double(Q(b, a));   // 10^-2
  assert(CGAL::is_finite(e));
  assert(e > 0.0099 && e < 0.0101);

  std::cout << "  large but representable: OK" << std::endl;
}
#endif // CGAL_USE_GMP

int main()
{
  std::cout << "Quotient<NT>::To_double (issues #1053, #1815)" << std::endl;

  // MP_Float has its own specialization of Real_embeddable_traits for
  // Quotient<MP_Float>, this is only a smoke test.
  std::cout << " Quotient<MP_Float>:" << std::endl;
  test_basic<CGAL::MP_Float>();

#ifdef CGAL_USE_GMP
  std::cout << " Quotient<Gmpz>:" << std::endl;
  test_basic<CGAL::Gmpz>();
  test_large_but_representable();
  test_out_of_range();
#else
  std::cout << " Skipping Gmpz tests, GMP is not available." << std::endl;
#endif

  std::cout << "All tests passed." << std::endl;
  return 0;
}
