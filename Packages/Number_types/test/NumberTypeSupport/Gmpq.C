#ifdef CGAL_USE_GMP

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <cassert>


typedef CGAL::Gmpz Gmpz;
typedef CGAL::Gmpq Gmpq;


 
int main(int, char **) {

  Gmpq q;
  Gmpq q1(12);
  Gmpq q2(3.1415);
  Gmpz z1(1), z2(2);
  Gmpq q3(z1);
  Gmpq q4(z1, z2);
  assert(q4.numerator() == Gmpz(1));
  assert(q4.denominator() == Gmpz(2));
  Gmpq q5(q4);

  Gmpq q6((signed long)1, 2);
  Gmpq q7((unsigned long)1, 2);
  return 0;
}
#else 
int main()
{
  return 0;
}
#endif //CGAL_USE_GMP
