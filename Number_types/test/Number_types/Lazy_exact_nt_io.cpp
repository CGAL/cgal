// Regression test for issue #135:
// Lazy_exact_nt operator<< used to_double(), losing precision and breaking
// round-trip save/load of exact kernel coordinates.
//
// The default operator<< still writes to_double(). CGAL::IO::set_exact_mode(os)
// switches it to write the exact value. For a rational backend (Exact_rational)
// that round-trips losslessly through operator>>; for CORE::Expr, whose operator<<
// writes a decimal approximation, exact mode is not lossless (see Test 8).

#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Exact_rational.h>

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_Expr.h>
#endif

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

typedef CGAL::Lazy_exact_nt<CGAL::Exact_rational> Lazy_nt;

int main()
{
  // Test 1: exact rational round-trip
  {
    std::cout << "Test 1: round-trip of 1/3" << std::endl;
    Lazy_nt a(CGAL::Exact_rational(1, 3));
    std::ostringstream oss;
    CGAL::IO::set_exact_mode(oss);
    oss << a;
    // Should output "1/3", not "0.333333..."
    std::istringstream iss(oss.str());
    Lazy_nt b;
    iss >> b;
    assert(iss);
    assert(a == b);
    std::cout << "  OK: wrote \"" << oss.str() << "\", read back equal" << std::endl;
  }

  // Test 2: integer value
  {
    std::cout << "Test 2: round-trip of 42" << std::endl;
    Lazy_nt a(42);
    std::ostringstream oss;
    CGAL::IO::set_exact_mode(oss);
    oss << a;
    std::istringstream iss(oss.str());
    Lazy_nt b;
    iss >> b;
    assert(iss);
    assert(a == b);
    std::cout << "  OK: wrote \"" << oss.str() << "\"" << std::endl;
  }

  // Test 3: negative rational
  {
    std::cout << "Test 3: round-trip of -5/7" << std::endl;
    Lazy_nt a(CGAL::Exact_rational(-5, 7));
    std::ostringstream oss;
    CGAL::IO::set_exact_mode(oss);
    oss << a;
    std::istringstream iss(oss.str());
    Lazy_nt b;
    iss >> b;
    assert(iss);
    assert(a == b);
    std::cout << "  OK: wrote \"" << oss.str() << "\"" << std::endl;
  }

  // Test 4: value that cannot be represented exactly as double
  // 1/3 + 1/7 = 10/21, to_double() would give 0.476190476190...
  {
    std::cout << "Test 4: round-trip of computed 1/3 + 1/7 = 10/21" << std::endl;
    Lazy_nt a(CGAL::Exact_rational(1, 3));
    Lazy_nt b(CGAL::Exact_rational(1, 7));
    Lazy_nt c = a + b;
    std::ostringstream oss;
    CGAL::IO::set_exact_mode(oss);
    oss << c;
    std::istringstream iss(oss.str());
    Lazy_nt d;
    iss >> d;
    assert(iss);
    assert(c == d);
    std::cout << "  OK: wrote \"" << oss.str() << "\"" << std::endl;
  }

  // Test 5: very large rational that would lose precision as double
  {
    std::cout << "Test 5: round-trip of large rational" << std::endl;
    CGAL::Exact_rational big;
    std::istringstream bigin("99999999999999999999/100000000000000000007");
    bigin >> big;
    assert(bigin);
    Lazy_nt a(big);
    std::ostringstream oss;
    CGAL::IO::set_exact_mode(oss);
    oss << a;
    std::istringstream iss(oss.str());
    Lazy_nt b;
    iss >> b;
    assert(iss);
    assert(a == b);
    std::cout << "  OK: wrote \"" << oss.str() << "\"" << std::endl;
  }

  // Test 6: zero
  {
    std::cout << "Test 6: round-trip of 0" << std::endl;
    Lazy_nt a(0);
    std::ostringstream oss;
    CGAL::IO::set_exact_mode(oss);
    oss << a;
    std::istringstream iss(oss.str());
    Lazy_nt b;
    iss >> b;
    assert(iss);
    assert(a == b);
    std::cout << "  OK: wrote \"" << oss.str() << "\"" << std::endl;
  }

  // Test 7: the default (no set_exact_mode) still writes a double, so existing
  // behaviour and user code are unchanged. Issue #135 must remain opt-in.
  {
    std::cout << "Test 7: default output is still to_double" << std::endl;
    Lazy_nt a(CGAL::Exact_rational(1, 3));
    std::ostringstream oss;
    oss << a; // no set_exact_mode: lossy, historical behaviour
    const std::string s = oss.str();
    assert(s.find('/') == std::string::npos); // a double, not "1/3"
    std::cout << "  OK: default wrote \"" << s << "\" (double, unchanged)" << std::endl;
  }

#ifdef CGAL_USE_CORE
  // Test 8: CORE::Expr regression. Exact mode is lossless only for backends whose
  // operator<< writes an exact representation. CORE::Expr::operator<< writes a
  // decimal approximation, so exact mode is NOT lossless for it. We check the
  // documented behaviour: the default still writes a double, and exact mode writes
  // exact()'s (approximate) representation, which parses.
  {
    typedef CGAL::Lazy_exact_nt<CORE::Expr> Lazy_core;
    std::cout << "Test 8: CORE::Expr exact mode is approximate, not lossless" << std::endl;
    Lazy_core a = Lazy_core(1) / Lazy_core(3); // 1/3

    std::ostringstream def;
    def << a; // default: a double, unchanged
    assert(def.str().find('/') == std::string::npos);

    std::ostringstream oss;
    CGAL::IO::set_exact_mode(oss);
    oss << a; // exact mode: CORE::Expr's decimal, parses but is not exact
    std::istringstream iss(oss.str());
    Lazy_core b;
    iss >> b;
    assert(iss);
    std::cout << "  OK: default \"" << def.str() << "\", exact-mode \""
              << oss.str() << "\" (approximate for CORE::Expr)" << std::endl;
  }
#endif

  std::cout << "All tests passed." << std::endl;
  return 0;
}

