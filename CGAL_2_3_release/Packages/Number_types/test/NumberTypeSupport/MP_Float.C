// Test program for the MP_Float class.
// Sylvain Pion.

#include <CGAL/basic.h>
#include <iostream>
#include <cassert>
#include <CGAL/MP_Float.h>
#include <CGAL/Random.h>
#include <CGAL/Quotient.h>

typedef CGAL::MP_Float MPF;

void test_equality(int i)
{
  assert(MPF(i) == MPF(i));
  assert(MPF(i) <= MPF(i));
  assert(MPF(i) >= MPF(i));
  assert(! (MPF(i) < MPF(i)));
  assert(! (MPF(i) > MPF(i)));
}

void bench(int loops)
{
  MPF j (13579);
  MPF k (9753);
  MPF l (-149);
  for(; loops>0; loops--)
  {
    MPF m (loops);
    MPF n = (j*k-l*m)*(j*m-k*l)*(j*j);
    (void) (n>MPF(0));
  }
}

MPF factoriel (short i)
{
  MPF r(1);
  for (short j=1; j<=i; j++)
    r = r * MPF(j);
  return r;
}

void print_test()
{
  for (int i=-2; i<3; i++)
    std::cout << MPF(i) << "   ";
  std::cout << std::endl;
  for (int j=-2; j<3; j++)
    std::cout << MPF(j+65536) << "      ";
  std::cout << std::endl;
}

int main(int argc, char **argv)
{
  CGAL::Quotient<CGAL::MP_Float> q1(1), q2(2);
  assert(q1+q1 == q2);

  int loops = argc > 1 ? atoi(argv[1]) : 100;
  bench(loops);

  std::cout.precision(20);

  for (; loops >= 0; loops--) {
    double d = CGAL::default_random.get_double();
    int exp = int((CGAL::default_random.get_double()-.5)*1024);
    d = ::ldexp(d, exp);
    // d = 7.34766e-140; // Crash encore sur PC...
    // d = 1.9696110926449043849e+124;
    // d = 7.3074228557478900057e+47;
    // std::cout << d << std::endl;
    // std::cout << MPF(d) << std::endl;
    // std::cout << CGAL::to_double(MPF(d)) << std::endl;
    if (CGAL::to_double(MPF(d)) != d)
      std::cerr << "CONVERSION ERROR with double : " << d << std::endl;
    // MPF z(d);
  }

  MPF a(0);
  int j=0;
  for (int i=0; i<10000; i++)
  {
#if 0
    std::cout << i << std::endl;
    std::cout << a << std::endl;
    std::cout << "exp = " << a.exp << "    v[" << a.v.size()<< "] = ";
    if (!a.is_zero())
      for (unsigned int j=0; j<a.v.size(); j++)
        std::cout << " " << a.v[j];
    std::cout << std::endl;
#endif
    a = a + MPF(i);
    j += i;
    CGAL_assertion(a == MPF(j));
    // std::cout << "a=" << a << std::endl << "MPF(j)=" << MPF(j) << std::endl;
  }

  std::cout << a << std::endl;
  CGAL_assertion ( a == MPF(5000*9999));

  MPF bb = factoriel(100);
  std::cout << "100! = " << bb << std::endl;
  
  MPF b = factoriel(10);
  std::cout << "10! = " << b << " =? 3628800 " << " =? " << CGAL::to_double(b);
  std::cout << std::endl;
  CGAL_assertion ( b == MPF(3628800));

  test_equality(0);
  test_equality(10);
  test_equality(-10);

  MPF c = -10*65536+30; // Note : limb==short is hardcoded here
  MPF d = c * 65536;

  CGAL_assertion(c < (c+1));
  CGAL_assertion(c > (c-1));
  CGAL_assertion(d < (d+1));
  CGAL_assertion(d > (d-1));

  MPF e =  7*65536 + 10;
  MPF f = 10*65536 + 7;
  CGAL_assertion(e < f);

  CGAL_assertion(MPF(-2) == -MPF(2));

  print_test();

  return 0;
}
