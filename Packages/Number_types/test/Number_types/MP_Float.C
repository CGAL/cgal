// Test program for the MP_Float class.
// Sylvain Pion.

#include <CGAL/basic.h>
#include <iostream>
#include <cassert>
#include <CGAL/MP_Float.h>
#include <CGAL/Random.h>
#include <CGAL/Quotient.h>

typedef CGAL::MP_Float       MPF;
typedef CGAL::Quotient<MPF>  QMPF;

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

void square_test()
{
  for (int i = 0; i<1000; ++i)
  {
    double d = CGAL::default_random.get_double();
    MPF D(d);
    CGAL_assertion(D*D == CGAL::square(D));
  }

  // Test case by Nico Kruithof.
  CGAL::MP_Float f = 32767.5;
  CGAL::MP_Float g = 32767.5;
  CGAL::MP_Float sqr = CGAL::square(f);
  assert(sqr == f*g);
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

void test_overflow_to_double()
{
  std::cout << "Tests if to_double(Quotient<MPF>) overflows or not."
            << std::endl;

  QMPF val = MPF(1)/2;
  for (int i=0; i<3000; ++i) {
    // std::cout << CGAL::to_double(val) << std::endl;
    // std::cout << val.numerator() << " , " << val.denominator() << std::endl;
    val = val * (1<<16);
    val = val / (1<<16);
  }
  assert(CGAL::to_double(val) == 0.5);
}

void test_overflow_to_interval()
{
  std::cout << "Tests if to_interval(Quotient<MPF>) overflows or not."
            << std::endl;

  QMPF val = MPF(1)/2;
  for (int i=0; i<3000; ++i) {
    // std::cout << CGAL::to_double(val) << std::endl;
    // std::cout << val.numerator() << " , " << val.denominator() << std::endl;
    val = val * (1<<16);
    val = val / (1<<16);
  }
  assert(CGAL::to_interval(val) == std::make_pair(0.5, 0.5));
}

void test_overflow_exponent()
{
  std::cout << "Testing overflow in exponent." << std::endl;
  QMPF a = 2;
  a = a/a;
  double loops = 0;
  while (loops < 200) // This should be able to loop forever.
  {
    a = a*a;
    assert((a+1)-1 == a);
    loops += 1;
    // std::cout << "loops = " << loops << std::endl;
  };
}

int main(int argc, char **argv)
{
  std::cout << "sizeof(long) = " << sizeof(long) << std::endl;
  std::cout << "sizeof(int) = " << sizeof(int) << std::endl;
  std::cout << "sizeof(short) = " << sizeof(short) << std::endl;
  std::cout << "sizeof(char) = " << sizeof(char) << std::endl;
  MPF z = int(65536);
  assert(z != 0);

  QMPF q1(1), q2(2);
  assert(q1+q1 == q2);

  int loops = argc > 1 ? CGAL_CLIB_STD::atoi(argv[1]) : 100;
  bench(loops);

  std::cout.precision(20);

  for (; loops >= 0; loops--) {
    double d = CGAL::default_random.get_double();
    int exp = int((CGAL::default_random.get_double()-.5)*1024);
    d = CGAL_CLIB_STD::ldexp(d, exp);
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
    std::cout << j << std::endl;
    std::cout << a << std::endl;
    std::cout << "exp = " << a.exp << "    v[" << a.v.size()<< "] = ";
    if (!a.is_zero())
      for (unsigned int j=0; j<a.v.size(); j++)
        std::cout << " " << a.v[j];
    std::cout << std::endl;
#endif
    a = a + MPF(i);
    j += i;
    assert(i == MPF(i));
#if 0
    std::cout << " array = "; print(std::cout, a); std::cout << std::endl;
    std::cout << i << std::endl;
    std::cout << j << std::endl;
    std::cout << a << std::endl;
    std::cout << "exp = " << a.exp << "    v[" << a.v.size()<< "] = ";
    if (!a.is_zero())
      for (unsigned int j=0; j<a.v.size(); j++)
        std::cout << " " << a.v[j];
    std::cout << std::endl;
#endif
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

  square_test();

  test_overflow_to_double();

  test_overflow_to_interval();

  test_overflow_exponent();

  return 0;
}
