// Regression test for issue #135:
// Lazy_exact_nt operator<< used to_double(), losing precision and breaking
// round-trip save/load of exact kernel coordinates.

#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Exact_rational.h>

#include <cassert>
#include <iostream>
#include <sstream>

typedef CGAL::Lazy_exact_nt<CGAL::Exact_rational> Lazy_nt;

int main()
{
  // Test 1: exact rational round-trip
  {
    std::cout << "Test 1: round-trip of 1/3" << std::endl;
    Lazy_nt a(CGAL::Exact_rational(1, 3));
    std::ostringstream oss;
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
    oss << a;
    std::istringstream iss(oss.str());
    Lazy_nt b;
    iss >> b;
    assert(iss);
    assert(a == b);
    std::cout << "  OK: wrote \"" << oss.str() << "\"" << std::endl;
  }

  std::cout << "All tests passed." << std::endl;
  return 0;
}

