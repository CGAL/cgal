// Regression test for https://github.com/CGAL/cgal/issues/8966
// and the related fixes for #1053 and #1815.
//
// Issue #8966: to_double(Quotient<Gmpzf>) was not monotonic — it was
// possible that a > b yet to_double(a) < to_double(b) even when both
// values lie in [0, double_max].
//
// Issues #1053 / #1815: the old To_double performed division in the
// ring type NT (not guaranteed for ring types) and mis-detected
// overflow by testing is_finite() on the NT rather than on the double.
//
// The rewritten To_double uses a fast-path (to_double(num)/to_double(den))
// when both fit in a finite double, and a slow-path (interval arithmetic
// on num and den, then midpoint) otherwise.  This test verifies
// correctness and monotonicity.

#include <cassert>
#include <iostream>

#include <CGAL/Quotient.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif

#include <CGAL/MP_Float.h>

// ---------- helpers ----------

template <typename QT>
void assert_monotonic(const QT& a, const QT& b,
                      const char* label)
{
  double da = CGAL::to_double(a);
  double db = CGAL::to_double(b);
  if (a > b) {
    if (!(da >= db)) {
      std::cerr << "MONOTONICITY FAILURE (" << label << "): "
                << "a > b but to_double(a) < to_double(b)" << std::endl;
      std::cerr << "  to_double(a)=" << da << "  to_double(b)=" << db
                << std::endl;
      assert(false);
    }
  } else if (a < b) {
    if (!(da <= db)) {
      std::cerr << "MONOTONICITY FAILURE (" << label << "): "
                << "a < b but to_double(a) > to_double(b)" << std::endl;
      std::cerr << "  to_double(a)=" << da << "  to_double(b)=" << db
                << std::endl;
      assert(false);
    }
  }
}

// ---------- #8966 exact reproducer ----------

#ifdef CGAL_USE_GMP
void test_issue_8966_reproducer()
{
  typedef CGAL::Quotient<CGAL::Gmpz> Q;

  std::cout << "  Testing #8966 exact reproducer..." << std::endl;

  // From the issue, the original Quotient<Gmpzf> values have the form
  //   m * 2^e
  // We reconstruct the same ratios as Quotient<Gmpz> by absorbing
  // the 2^e factors into numerator and denominator:
  //
  //   a = (num_a * 2^-379) / (den_a * 2^-299)
  //     = num_a / (den_a * 2^80)
  //
  //   b = (num_b * 2^-215) / (den_b * 2^-141)
  //     = num_b / (den_b * 2^74)

  CGAL::Gmpz num_a(
    "495465331884540240762104420639278860096116250934219950844827973437034498625");
  CGAL::Gmpz den_a(
    "2245694908428994193174821578352066558595985337089");
  CGAL::Gmpz num_b(
    "3447327498334006902041169");
  CGAL::Gmpz den_b("1");

  CGAL::Gmpz two_80("1208925819614629174706176");   // 2^80
  CGAL::Gmpz two_74("18889465931478580854784");      // 2^74

  Q a(num_a, den_a * two_80);
  Q b(num_b, den_b * two_74);

  double da = CGAL::to_double(a);
  double db = CGAL::to_double(b);
  std::cout << "    a > b: " << (a > b ? "true" : "false") << std::endl;
  std::cout << "    to_double(a) = " << da << std::endl;
  std::cout << "    to_double(b) = " << db << std::endl;

  assert(a > b);
  assert_monotonic(a, b, "issue-8966-reproducer");
}
#endif // CGAL_USE_GMP

// ---------- generic Quotient<NT> tests ----------

template <typename NT>
void test_to_double_basic()
{
  typedef CGAL::Quotient<NT> Q;

  // Zero.
  {
    Q q(NT(0), NT(1));
    assert(CGAL::to_double(q) == 0.0);
  }

  // Denominator is 1.
  {
    Q q(NT(42), NT(1));
    assert(CGAL::to_double(q) == 42.0);
  }

  // Simple fraction.
  {
    Q q(NT(1), NT(2));
    double d = CGAL::to_double(q);
    assert(d == 0.5);
  }

  // Negative.
  {
    Q q(NT(-7), NT(3));
    double d = CGAL::to_double(q);
    assert(d < 0.0);
    // Should be close to -2.333...
    assert(d > -2.4 && d < -2.3);
  }
}

template <typename NT>
void test_to_double_monotonicity()
{
  typedef CGAL::Quotient<NT> Q;

  // Construct pairs where a > b and check to_double(a) >= to_double(b).

  // Pair 1: values near each other.
  {
    Q a(NT(100), NT(3));   // 33.333...
    Q b(NT(100), NT(4));   // 25.0
    assert(a > b);
    assert_monotonic(a, b, "near-values");
  }

  // Pair 2: both components large.
  {
    NT big1(1);
    NT big2(1);
    for (int i = 0; i < 40; ++i) {
      big1 = big1 * NT(10);
      big2 = big2 * NT(10);
    }
    big1 = big1 + NT(1);  // slightly larger

    Q a(big1, NT(1));
    Q b(big2, NT(1));
    assert(a > b);
    assert_monotonic(a, b, "large-values");
  }

  // Pair 3: overflow regime — both to_double(num) and to_double(den) overflow.
  // 10^400 > DBL_MAX ≈ 1.8e308, so to_double() returns ±inf and the
  // slow path (interval arithmetic) is exercised.
  //
  // Note: for generic ring types without a To_double specialization
  // (none of the standard CGAL number types fall into this category),
  // the slow path returns a very imprecise result when both num and den
  // overflow — but monotonicity is still preserved.
  {
    NT huge(1);
    for (int i = 0; i < 400; ++i)
      huge = huge * NT(10);

    Q a(huge + NT(1), huge);          // slightly > 1
    Q b(huge, huge + NT(1));          // slightly < 1
    assert(a > b);
    assert_monotonic(a, b, "overflow-regime");
  }
}

int main()
{
  std::cout.precision(20);
  std::cout << "Testing Quotient<NT>::To_double (issues #1053, #1815, #8966)..."
            << std::endl;

  // Test with MP_Float (ring type — no operator/).
  // Note: Real_embeddable_traits<Quotient<MP_Float>> has a specialization
  // in MP_Float.h that bypasses the generic To_double in Quotient.h.
  // This section tests that specialization, not the generic code path.
  std::cout << "  Quotient<MP_Float>:" << std::endl;
  test_to_double_basic<CGAL::MP_Float>();
  test_to_double_monotonicity<CGAL::MP_Float>();

#ifdef CGAL_USE_GMP
  // Test with Gmpz.
  // Real_embeddable_traits<Quotient<Gmpz>> uses the generic To_double
  // in Quotient.h (the Algebraic_structure_traits specialization in
  // Gmpz.h does not affect CGAL::to_double()).
  std::cout << "  Quotient<Gmpz>:" << std::endl;
  test_to_double_basic<CGAL::Gmpz>();
  test_to_double_monotonicity<CGAL::Gmpz>();

  // Issue #8966 exact reproducer.
  test_issue_8966_reproducer();
#else
  std::cout << "  Skipping Gmpz / Gmpzf tests (GMP not available)." << std::endl;
#endif

  std::cout << "All tests passed." << std::endl;
  return 0;
}
