#ifdef CGAL_USE_GMP

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <cassert>


typedef CGAL::Gmpz Gmpz;
typedef CGAL::Gmpq Gmpq;


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
  return 0;
}
#else 
int main()
{
  return 0;
}
#endif //CGAL_USE_GMP
