#ifdef CGAL_USE_GMP

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>

#ifdef CGAL_USE_GMPXX
#  include <CGAL/gmpxx.h>
#endif

#include <CGAL/Interval_nt.h>
#include <cassert>


typedef CGAL::Gmpz Gmpz;
typedef CGAL::Gmpq Gmpq;

template < typename NT >
void test_overflow_to_interval(const NT&)
{
  std::cout << "Tests if to_interval() handles overflow nicely or not."
            << std::endl;

  NT val = 2;
  for (int i=0; i<20; ++i) {
    val = val*val;
    CGAL::Interval_nt<> inter = CGAL::to_interval(val);
    CGAL::Interval_nt<> minter = CGAL::to_interval(-val);
    assert(CGAL::is_valid(inter));
    assert(CGAL::is_finite(inter.inf()));
    assert(CGAL::is_valid(minter));
    assert(CGAL::is_finite(minter.sup()));
  }
}

void test_overflow_to_double()
{
  std::cout << "Tests if to_double(Gmpq) overflows or not." << std::endl;

  Gmpq val = Gmpq(1)/2;
  for (int i=0; i<3000; ++i) {
    // std::cout << CGAL::to_double(val) << std::endl;
    // std::cout << val.numerator() << " , " << val.denominator() << std::endl;
    val = val * (1<<16);
    val = val / (1<<16);
  }
  assert(CGAL::to_double(val) == 0.5);
}

int main() {

  Gmpq q;
  Gmpq q1(12);
  Gmpq q2(3.1415);
  Gmpz z1(1), z2(2);
  Gmpq q3(z1);
  Gmpq q4(z1, z2);
  assert(q4.numerator() == Gmpz(1));
  assert(q4.denominator() == Gmpz(2));
  Gmpq q5(q4);

  Gmpq qi1(0,1); // int int
  Gmpq qi2(3, -3);
  assert(qi2.numerator() == Gmpz(-1)); // because Gmpq is canonicalized
  Gmpq qi3(-3, 3);
  assert(qi2 == qi3);

  Gmpq q6((signed long)1, (unsigned long)2);
  Gmpq q7((unsigned long)1, (unsigned long)2);

  Gmpq qs("-100/-111", 2);
  assert(qs.numerator() == Gmpz(4));
  assert(qs.denominator() == Gmpz(7));

  Gmpq qs2("100/1");
  assert(qs2.numerator() == Gmpz(100));
  assert(qs2.denominator() == Gmpz(1));

  test_overflow_to_double();
  test_overflow_to_interval(Gmpz());
  test_overflow_to_interval(Gmpq());
#ifdef CGAL_USE_GMPXX
  test_overflow_to_interval(mpz_class());
  test_overflow_to_interval(mpq_class());
#endif

  return 0;
}
#else 
int main()
{
  return 0;
}
#endif //CGAL_USE_GMP
