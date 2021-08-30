// Test program for simplest_rational_in_interval().

#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <cassert>

#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#  include <CGAL/Gmpq.h>
#endif

#ifdef CGAL_USE_GMPXX
#  include <CGAL/gmpxx.h>
#endif

#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_rational.h>
#endif

#include <CGAL/simplest_rational_in_interval.h>
#include <CGAL/to_rational.h>

template <class Q>
void test_it()
{
  Q q = CGAL::simplest_rational_in_interval<Q>(-0.1, 0.1);
  assert(CGAL_NTS is_zero(q));

  double l = 3.1415, h = 3.1416;
  q = CGAL::simplest_rational_in_interval<Q>(l, h);
  assert(l <= CGAL_NTS to_double(q));
  assert(CGAL_NTS to_double(q) <= h);

  double d = 1234.56789;
  q = CGAL:: to_rational<Q>(d);
  assert(CGAL_NTS to_double(q) == d);
}

int main() {
  test_it<CGAL::Quotient<CGAL::MP_Float> >();

#ifdef CGAL_USE_GMP
  test_it<CGAL::Quotient<CGAL::Gmpz> >();
  test_it<CGAL::Gmpq>();
#endif

#ifdef CGAL_USE_GMPXX
  test_it<mpq_class>();
#endif

#ifdef CGAL_USE_LEDA
  test_it<leda_rational>();
#endif

  return 0;
}
