

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Quotient.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
//#include <CGAL/gmpxx.h>
#endif // CGAL_USE_GMP

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_rational.h>
#endif // CGAL_USE_LEDA


#include <CGAL/simplest_rational_in_interval.h>
#include <CGAL/to_rational.h>


#ifdef CGAL_USE_GMP
typedef CGAL::Gmpz Gmpz;
typedef CGAL::Gmpq Gmpq;
typedef CGAL::Quotient<CGAL::Gmpz> Quot;
#else
typedef CGAL::Quotient<int> Quot;
#endif

template <class Q>
void test_it()
{
  Q q = CGAL::simplest_rational_in_interval<Q>(-0.1, 0.1);
  assert(CGAL::NTS::is_zero(q));

  double l = 3.1415, h = 3.1416;
  q = CGAL::simplest_rational_in_interval<Q>(l, h);
  assert(l <= CGAL::to_double(q));
  assert(CGAL::to_double(q) <= h);
  
  double d = 1234.56789;
  q = CGAL::to_rational<Q>(d);
  assert(CGAL::to_double(q) == d);
}
 
int main(int, char **) {
#ifdef CGAL_USE_GMP
  test_it<Gmpq>();
#endif //CGAL_USE_GMP

  //#ifdef CGAL_USE_GMPXX
  //test_it<mpq_class>();
  //#endif

  test_it<Quot>();


#ifdef CGAL_USE_LEDA
  test_it<leda_rational>();
#endif

  return 0;
}


