// Regression test for https://github.com/CGAL/cgal/issues/5990
//
// to_interval(cpp_int) threw "Cannot convert a non-finite number to an
// integer" for values whose double conversion overflows to infinity.
// The old implementation called x.compare(double), which Boost.MP
// implements via convert_to<number>, failing on ±infinity.
//
// The fix (Boost_MP_internal::to_interval) uses bit-level extraction
// instead.  This test ensures the fix does not regress.

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

#include <CGAL/config.h>

#ifdef CGAL_USE_BOOST_MP
#include <CGAL/boost_mp.h>
#endif

int main()
{
#ifdef CGAL_USE_BOOST_MP
  typedef boost::multiprecision::cpp_int I;
  const double inf = std::numeric_limits<double>::infinity();

  std::cout << "Testing to_interval(cpp_int) for large numbers (issue #5990)..."
            << std::endl;

  // --- Exact reproducer from the issue ---
  // i = 2^10000, which overflows double to +inf.
  // The old code threw std::runtime_error here.
  {
    I i = 1;
    i <<= 10000;
    auto interval = CGAL::to_interval(i);
    double lo = interval.first;
    double hi = interval.second;
    std::cout << "  2^10000: [" << lo << ", " << hi << "]" << std::endl;
    // The value overflows double, so the tight interval is [DBL_MAX, +inf].
    assert(lo == (std::numeric_limits<double>::max)());
    assert(hi == inf);
  }

  // --- Negative large value ---
  {
    I i = -1;
    i <<= 10000;
    auto interval = CGAL::to_interval(i);
    double lo = interval.first;
    double hi = interval.second;
    std::cout << "  -2^10000: [" << lo << ", " << hi << "]" << std::endl;
    // Symmetric to the positive case: the tight interval is [-inf, -DBL_MAX].
    assert(lo == -inf);
    assert(hi == std::numeric_limits<double>::lowest());
  }

  // --- Moderately large value: 2^1024 (just above double_max ≈ 2^1023.999) ---
  {
    I i = 1;
    i <<= 1024;
    auto interval = CGAL::to_interval(i);
    double lo = interval.first;
    double hi = interval.second;
    std::cout << "  2^1024: [" << lo << ", " << hi << "]" << std::endl;
    // 2^1024 > DBL_MAX, so it also overflows: tight interval is [DBL_MAX, +inf].
    assert(lo == (std::numeric_limits<double>::max)());
    assert(hi == inf);
  }

  // --- Value that fits in double: 2^52 ---
  {
    I i = 1;
    i <<= 52;
    auto interval = CGAL::to_interval(i);
    double lo = interval.first;
    double hi = interval.second;
    double expected = std::ldexp(1.0, 52);
    std::cout << "  2^52: [" << lo << ", " << hi << "]" << std::endl;
    assert(lo == expected);
    assert(hi == expected);
  }

  // --- Zero ---
  {
    I i = 0;
    auto interval = CGAL::to_interval(i);
    std::cout << "  0: [" << interval.first << ", " << interval.second << "]"
              << std::endl;
    assert(interval.first == 0.0);
    assert(interval.second == 0.0);
  }

  // --- Large value that is exactly representable: 2^100 ---
  {
    I i = 1;
    i <<= 100;
    auto interval = CGAL::to_interval(i);
    double lo = interval.first;
    double hi = interval.second;
    double expected = std::ldexp(1.0, 100);
    std::cout << "  2^100: [" << lo << ", " << hi << "]" << std::endl;
    assert(lo == expected);
    assert(hi == expected);
  }

  std::cout << "All tests passed." << std::endl;

#else
  std::cout << "Test skipped (CGAL_USE_BOOST_MP not defined)." << std::endl;
#endif

  return 0;
}
