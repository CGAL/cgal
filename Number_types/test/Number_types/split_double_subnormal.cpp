// Regression test for https://github.com/CGAL/cgal/issues/5982
//
// Building CGAL::Quotient<Integer> from a subnormal (denormal) double used to
// abort: Split_double relied on split_numerator_denominator(), whose
// denominator overflows to infinity for subnormals (denorm_min needs a
// denominator of 2^1074), so CGAL_postcondition(d == num/den) failed. The fix
// decomposes the double with split_mantissa_exponent() instead. This test
// builds Quotient<Integer> from several subnormal doubles and checks the
// result is exact, for every available big-integer backend.

#include <CGAL/config.h>
#include <CGAL/Quotient.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif
#ifdef CGAL_USE_GMPXX
#include <CGAL/mpz_class.h>
#endif
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#endif
#ifdef CGAL_USE_BOOST_MP
#include <CGAL/boost_mp.h>
#endif

// x *= 2^k (k >= 0), using only addition so it works for every integer type
// (including leda_integer, which does not expose operator<<=).
template <typename NT>
void mul_pow2(NT& x, int k)
{
  for (int i = 0; i < k; ++i) x = x + x;
}

// Exact reference decomposition of d into num/den, independent of the
// library's Split_double. The integer is built from the double mantissa
// (|mantissa| < 2^53, hence exact), as every Split_double specialization does.
template <typename NT>
CGAL::Quotient<NT> reference_quotient(double d)
{
  if (d == 0.0) return CGAL::Quotient<NT>(NT(0), NT(1));
  int exp;
  const double frac = std::frexp(d, &exp);      // d == frac * 2^exp
  const double mantissa = std::ldexp(frac, 53); // exact signed integer
  int e = exp - 53;                             // d == mantissa * 2^e
  NT num(mantissa), den(1);
  if (e >= 0) mul_pow2(num, e);
  else        mul_pow2(den, -e);
  return CGAL::Quotient<NT>(num, den);
}

template <typename NT>
void test(const char* name)
{
  typedef std::numeric_limits<double> lim;
  const double values[] = {
    lim::denorm_min(),                 // 2^-1074, the issue #5982 reproducer
    -lim::denorm_min(),
    3.0 * lim::denorm_min(),
    std::nextafter(lim::min(), 0.0),   // largest subnormal
    lim::min() / 2.0,                  // a subnormal
    lim::min(),                        // smallest normal (boundary)
    0.0, 1.0, -1.0, 0.5, -0.25, 42.0
  };
  for (double d : values)
  {
    CGAL::Quotient<NT> q(d);           // aborted here for subnormals before the fix
    assert(q == reference_quotient<NT>(d));
  }

  // Independent anchor that does not go through frexp: denorm_min * 2^1074 == 1.
  {
    CGAL::Quotient<NT> q(lim::denorm_min());
    NT p(1); mul_pow2(p, 1074);
    assert(q * p == CGAL::Quotient<NT>(NT(1)));
  }

  std::cout << "  " << name << ": ok" << std::endl;
}

int main()
{
  std::cout << "Testing Split_double for subnormal doubles (issue #5982)..."
            << std::endl;
#ifdef CGAL_USE_GMP
  test<CGAL::Gmpz>("Gmpz");
#endif
#ifdef CGAL_USE_GMPXX
  test<mpz_class>("mpz_class");
#endif
#ifdef CGAL_USE_LEDA
  test<leda_integer>("leda_integer");
#endif
#ifdef CGAL_USE_BOOST_MP
  test<boost::multiprecision::cpp_int>("cpp_int");
#endif
  std::cout << "All tests passed." << std::endl;
  return 0;
}
