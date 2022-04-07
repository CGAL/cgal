// Test program for the MP_Float class.
// Sylvain Pion.

#include <CGAL/config.h>
#include <CGAL/exceptions.h>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <CGAL/MP_Float.h>
#include <CGAL/Random.h>
#include <CGAL/Quotient.h>
#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpq.h>
#endif

typedef CGAL::MP_Float       MPF;
typedef CGAL::Quotient<MPF>  QMPF;

double non_zero_double(){
 double d;
 do {
  d = CGAL::get_default_random().get_double();
  if(d ==0) {
          std::cout << "generated zero" << std::endl;
  }
 }while(d==0);
 return d;
}

void test_is_integer()
{
  std::cout << "Testing is_integer()" << std::endl;
  assert(! is_integer(MPF(0.5)));
  assert(! is_integer(MPF(0.25)));
  assert(! is_integer(MPF(0.1)));
  assert(! is_integer(MPF(1e-100)));
  assert(! is_integer(MPF(1.5)));
  assert(! is_integer(MPF(15.1)));
  assert(! is_integer(MPF(1e100)+0.5));
  assert(! is_integer(MPF(-0.5)));
  assert(! is_integer(MPF(-0.25)));
  assert(! is_integer(MPF(-0.1)));
  assert(! is_integer(MPF(-1e-100)));
  assert(! is_integer(MPF(-1.5)));
  assert(! is_integer(MPF(-15.1)));
  assert(! is_integer(MPF(-1e100)+0.5));

  assert(is_integer(MPF(0)));
  assert(is_integer(MPF(1)));
  assert(is_integer(MPF(2)));
  assert(is_integer(MPF(1e100)));
  assert(is_integer(MPF(1e100)+1));
  assert(is_integer(MPF(-0)));
  assert(is_integer(MPF(-1)));
  assert(is_integer(MPF(-2)));
  assert(is_integer(MPF(-1e100)));
  assert(is_integer(MPF(-1e100)+1));
}

void test_integral_division()
{
  std::cout << "Testing integral_division()" << std::endl;

  // Let's pick 2 random values, multiply them, divide and check.

  MPF tmp = integral_division(MPF(0), MPF(1));

  for (int i = 0; i < 10000; ++i) {
    MPF d = non_zero_double();
    MPF n = non_zero_double();

    // We test with up to 5 chunks for the denominator.
    if ((i%5) < 1) d *= non_zero_double();
    if ((i%5) < 2) d *= non_zero_double();
    if ((i%5) < 3) d *= non_zero_double();
    if ((i%5) < 4) d *= non_zero_double();

    // We test with up to 3 chunks for the numerator.
    if ((i%3) < 1) n *= non_zero_double();
    if ((i%3) < 2) n *= non_zero_double();

    // Try to change the signs
    if ((i%7)  == 0) d = -d;
    if ((i%11) == 0) n = -n;

    // Scale both numerator and denominator to test overflow/underflow
    // situations.
    d.rescale(107 * ((i%19) - 10));
    n.rescale(102 * ((i%33) - 16));

    MPF nd = d * n;
    assert(nd == n * d);
    //std::cout << "d = " << d << std::endl;
    //std::cout << "n = " << n << std::endl;
    assert( integral_division(nd, n) == d);
    assert( integral_division(nd, d) == n);
    assert( CGAL::divides(n,nd) );
    assert( CGAL::divides(d,nd) );
  }

  assert( ! CGAL::divides(MPF(3), MPF(1)) );
  assert( ! CGAL::divides(MPF(7), MPF(2)) );
  // test if we're lucky :)
  assert( ! CGAL::divides(MPF(non_zero_double()), MPF(non_zero_double())) );
}

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
    double d = CGAL::get_default_random().get_double();
    MPF D(d);
    assert(D*D == CGAL_NTS square(D));
  }

  // Test case by Nico Kruithof.
  CGAL::MP_Float f = 32767.5;
  CGAL::MP_Float g = 32767.5;
  CGAL::MP_Float sqr = CGAL_NTS square(f);
  assert(sqr == f*g);

  // Another test case.
  CGAL::MP_Float M_16 = 1.0/(1<<16);
  CGAL::MP_Float M_32 = M_16*M_16;
  CGAL::MP_Float M_48 = M_32*M_16;
  CGAL::MP_Float M_64 = M_48*M_16;
  CGAL::MP_Float m = -32768 * M_64 -4795 * M_48 +10121 * M_32 -32768 * M_16;
  // square(m) is wrongly computed.
  // For now, I disable square(MP_Float).
  CGAL::MP_Float s = CGAL_NTS square(m);
  assert(s*4 == CGAL_NTS square(2*m));
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

  QMPF val = QMPF(1)/2;
  for (int i=0; i<3000; ++i) {
    // std::cout << CGAL::to_double(val) << std::endl;
    // std::cout << val.numerator() << " , " << val.denominator() << std::endl;
    val = val * (1<<16);
    val = val / (1<<16);
  }
  assert(CGAL_NTS to_double(val) == 0.5);
}

void test_overflow_to_interval()
{
  std::cout << "Tests if to_interval(Quotient<MPF>) overflows or not."
            << std::endl;

  QMPF val = QMPF(1)/2;
  for (int i=0; i<3000; ++i) {
    // std::cout << CGAL::to_double(val) << std::endl;
    // std::cout << val.numerator() << " , " << val.denominator() << std::endl;
    val = val * (1<<16);
    val = val / (1<<16);
  }
  assert(CGAL_NTS to_interval(val) == std::make_pair(0.5, 0.5));

  // Reported by Pedro to_interval(10^400) returns [+inf;+inf] :
  MPF m = 10;
  for (int j=0; j<400; ++j)
    m = 10 * m;
  CGAL::Interval_nt<> inter = CGAL_NTS to_interval(m);
  std::cout << inter << std::endl;
  assert(CGAL_NTS is_finite(inter.inf()));
  assert(!CGAL_NTS is_finite(inter.sup()));

  // Testing for underflow as well.
  MPF mm = 0.1;
  for (int k=0; k<1000; ++k)
    mm = mm * 0.125;
  CGAL::Interval_nt<> inter2 = CGAL_NTS to_interval(mm);
  std::cout << inter2 << std::endl;
  assert(inter2.inf() == 0);
  assert(inter2.sup() != 0);
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
  using CGAL::MP_Float;

  const unsigned        log_limb         = 8 * sizeof(MP_Float::limb);
  const MP_Float::limb2 base             = 1 << log_limb;
  const MP_Float::V::size_type limbs_per_double = 2 + 53/log_limb;

  const double trunc_max = double(base)*(base/2-1)/double(base-1);
  const double trunc_min = double(-base)*(base/2)/double(base-1);

  std::cout.precision(20);

  std::cout << " base = " << base << std::endl;
  std::cout << " trunc_min = " << trunc_min << std::endl;
  std::cout << " trunc_max = " << trunc_max << std::endl;
  std::cout << " limbs_per_double = " << limbs_per_double << std::endl;

  // Bug report by Stefan Uhrig on 18/06/07:
  double bug = 0.49999237058324297;
  std::cout << "bug = " << bug << std::endl;
  CGAL::MP_Float mpf(bug);

  std::cout << "sizeof(long) = " << sizeof(long) << std::endl;
  std::cout << "sizeof(int) = " << sizeof(int) << std::endl;
  std::cout << "sizeof(short) = " << sizeof(short) << std::endl;
  std::cout << "sizeof(char) = " << sizeof(char) << std::endl;
  MPF z = int(65536);
  assert(z != 0);

  QMPF q1(1), q2(2);
  assert(q1+q1 == q2);

  int loops = argc > 1 ? std::atoi(argv[1]) : 100;
  bench(loops);

  std::cout.precision(20);

  std::cout << "Checking MP_Float(float) constructor." << std::endl;
  for (int i = 0; i < loops; ++i) {
    float d = (float)CGAL::get_default_random().get_double();
    int exp = int((CGAL::get_default_random().get_double()-.5)*256);
    d = std::ldexp(d, exp);
    // std::cout << d << std::endl;
    // std::cout << MPF(d) << std::endl;
    // std::cout << CGAL_NTS to_double(MPF(d)) << std::endl;
    if (CGAL_NTS to_double(MPF(d)) != (double) d) {
      std::cerr << "CONVERSION ERROR with float : " << d << std::endl;
      std::abort();
    }
  }

  std::cout << "Checking MP_Float(double) constructor." << std::endl;
  MPF y = 0.5000000000000001; // see bug-report on cgal-discuss (2006-06-23).
  for (int i = 0; i < loops; ++i) {
    double d = CGAL::get_default_random().get_double();
    int exp = int((CGAL::get_default_random().get_double()-.5)*1024);
    d = std::ldexp(d, exp);
    // std::cout << d << std::endl;
    // std::cout << MPF(d) << std::endl;
    // std::cout << CGAL_NTS to_double(MPF(d)) << std::endl;
    if (CGAL_NTS to_double(MPF(d)) != d) {
      std::cerr << "CONVERSION ERROR with double : " << d << std::endl;
      std::abort();
    }
  }

  std::cout << "Checking MP_Float(long double) constructor." << std::endl;
  for (int i = 0; i < loops; ++i) {
    long double d = CGAL::get_default_random().get_double();
    d = d*d; // to get more bits
    int exp = int((CGAL::get_default_random().get_double()-.5)*1024);
    d = d * std::ldexp(1.0, exp);
    //std::cout << d << std::endl;
    //std::cout << MPF(d) << std::endl;
    //std::cout << CGAL_NTS to_double(MPF(d)) << std::endl;

    // Here we have to use a lesser test.
    MPF res = d;
    std::pair<double, double> ia = CGAL_NTS to_interval(res);
    if ( MPF(ia.first) > res || MPF(ia.second) < res) {
      std::cerr << "CONVERSION ERROR with long double : " << d << std::endl;
      std::abort();
    }
  }

  // Test-cases for specific bugs found :
  {
    double d = 0.4999974472643795;
    assert( d == CGAL_NTS to_double(MPF(d)) );
  }
  {
    double d = 7.3074228557478900057e+47;
    assert( d == CGAL_NTS to_double(MPF(d)) );
  }
  {
    double d = 1.9696110926449043849e+124;
    assert( d == CGAL_NTS to_double(MPF(d)) );
  }
  {
    double d = 7.34766e-140; // Crashed on PC...
    assert( d == CGAL_NTS to_double(MPF(d)) );
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
    assert(a == MPF(j));
    // std::cout << "a=" << a << std::endl << "MPF(j)=" << MPF(j) << std::endl;
  }

  std::cout << a << std::endl;
  assert ( a == MPF(5000*9999));

  MPF bb = factoriel(100);
  std::cout << "100! = " << bb << std::endl;

  MPF b = factoriel(10);
  std::cout << "10! = " << b << " =? 3628800 " << " =? " << CGAL_NTS to_double(b);
  std::cout << std::endl;
  assert ( b == MPF(3628800));

  test_equality(0);
  test_equality(10);
  test_equality(-10);

  MPF c = -10*65536+30; // Note : limb==short is hardcoded here
  MPF d = c * 65536;

  assert(c < (c+1));
  assert(c > (c-1));
  assert(d < (d+1));
  assert(d > (d-1));

  MPF e =  7*65536 + 10;
  MPF f = 10*65536 + 7;
  assert(e < f);

  assert(MPF(-2) == -MPF(2));

  print_test();

  square_test();

  test_overflow_to_double();

  test_overflow_to_interval();

  test_overflow_exponent();

#ifdef CGAL_USE_GMP
  double dt = -135.9682;
  MPF mdt = dt;
  CGAL::Gmpq q = mdt.to_rational<CGAL::Gmpq>();
  assert(q == CGAL::Gmpq(dt));
#endif // CGAL_USE_GMP

  test_integral_division();

  test_is_integer();

  return 0;
}
