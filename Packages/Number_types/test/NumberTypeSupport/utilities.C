#include <CGAL/basic.h>
#include <CGAL/Quotient.h> 
#include <CGAL/MP_Float.h> 
#include <CGAL/Lazy_exact_nt.h> 
#include <CGAL/Fixed_precision_nt.h> 

#ifndef CGAL_CFG_MATCHING_BUG_2
#include <CGAL/Filtered_exact.h> 
#endif // CGAL_CFG_MATCHING_BUG_2

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif // CGAL_USE_GMP

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>
#endif // CGAL_USE_LEDA

typedef CGAL::Quotient<CGAL::MP_Float>            QMPF;

// Exclude Filtered_exact tests with VC++ (lack of partial specialization)
#ifndef CGAL_CFG_MATCHING_BUG_2

typedef CGAL::Filtered_exact<double, QMPF>        FEDQ;
#ifdef CGAL_USE_GMP
typedef CGAL::Filtered_exact<int, CGAL::Gmpz>     FEIG;
typedef CGAL::Filtered_exact<double, CGAL::Gmpz>  FEDG;
#endif // CGAL_USE_GMP

#ifdef CGAL_USE_LEDA
typedef CGAL::Filtered_exact<double, leda_real>   FEDR;
#endif // CGAL_USE_LEDA

#endif // CGAL_CFG_MATCHING_BUG_2

template < class NT >
bool
test_sqrt(NT, CGAL::Tag_true)
{
  // sqrt
  NT sixteen(16);
  NT four(4);
  std::cerr << "  sqrt()" << std::endl;
  if (CGAL_NTS sqrt(sixteen) != four) return false;
  CGAL::Sqrt<NT> s;
  if (s(sixteen) != four) return false;

  return true;
}

template < class NT >
bool
test_sqrt(NT, CGAL::Tag_false)
{
  return true;
}

template < class NT >
bool
test_gcd(NT z, CGAL::Tag_true)
{
  // div
  NT eleven(11);
  NT three(3);
  std::cerr << "  div()" << std::endl;
  if (CGAL_NTS div(eleven, three) != three) return false;
  CGAL::Div<NT> d;
  if (d(eleven, three) != three) return false;
  
  // gcd
  NT x(2737);
  NT y(5083);
  std::cerr << "  gcd()" << std::endl;
  CGAL::Gcd<NT> gc;
  if (CGAL_NTS gcd(x, y) != NT(391)) return false;
  if (gc(x, y) != NT(391)) return false;

  typedef typename CGAL::Number_type_traits<NT>::Has_sqrt Has_sqrt;
  Has_sqrt has_sqrt = Has_sqrt();
  return test_sqrt(x, has_sqrt);
}

template < class NT >
bool
test_gcd(NT x, CGAL::Tag_false)
{
  typedef typename CGAL::Number_type_traits<NT>::Has_sqrt Has_sqrt;
  Has_sqrt has_sqrt = Has_sqrt();
  return test_sqrt(x, has_sqrt);
}

template < class NT >
bool
test_basic(NT x)
{
  NT zero(0);
  NT one(1);
  NT mone(-one);

  // is_zero
  std::cerr << "  is_zero()" << std::endl;
  if (!CGAL_NTS is_zero(zero)) return false;
  CGAL::Is_zero<NT> iz;
  if (!iz(zero)) return false;

  // is_one
  std::cerr << "  is_one()" << std::endl;
  if (!CGAL_NTS is_one(one)) return false;
  CGAL::Is_one<NT> io;
  if (!io(one)) return false;

  // is_positive
  std::cerr << "  is_positive()" << std::endl;
  if (!CGAL_NTS is_positive(one)) return false;
  CGAL::Is_positive<NT> ip;
  if (!ip(one)) return false;

  // is_negative
  std::cerr << "  is_negative()" << std::endl;
  if (CGAL_NTS is_negative(one)) return false;
  CGAL::Is_negative<NT> in;
  if (in(one)) return false;

  // sign
  std::cerr << "  sign()" << std::endl;
  CGAL::Sgn<NT> sg;
  if (CGAL_NTS sign(one) != CGAL::POSITIVE) return false;
  if (sg(one) != CGAL::POSITIVE) return false;
  if (CGAL_NTS sign(zero) != CGAL::ZERO) return false;
  if (sg(zero) != CGAL::ZERO) return false;

  // abs
  std::cerr << "  abs()" << std::endl;
  CGAL::Abs<NT> ab;
  if (CGAL_NTS abs(mone) != one &&
      CGAL_NTS abs(mone) != mone) // unsigned types :-)
    return false;
  if (ab(mone) != one && ab(mone) != mone) 
    return false;
  if (ab(mone) != CGAL_NTS abs(mone) || ab(one) != CGAL_NTS abs(one))
    return false;

  // compare
  std::cerr << "  compare()" << std::endl;
  CGAL::Compare<NT> co;
  if (CGAL_NTS compare(one, zero) != CGAL::LARGER) return false;
  if (CGAL_NTS compare(one, one) != CGAL::EQUAL) return false;
  if (CGAL_NTS compare(zero, one) != CGAL::SMALLER) return false;
  if (co(one, zero) != CGAL::LARGER) return false;
  if (co(one, one) != CGAL::EQUAL) return false;
  if (co(zero, one) != CGAL::SMALLER) return false;

  // square
  std::cerr << "  square()" << std::endl;
  NT two(2);
  NT four(4);
  CGAL::Square<NT> sq;
  if (CGAL_NTS square(two) != four) return false;
  if (sq(two) != four) return false;

  typedef typename CGAL::Number_type_traits<NT>::Has_gcd Has_gcd;
  Has_gcd has_gcd = Has_gcd();
  return test_gcd(x, has_gcd);
}

#define TESTIT(T,N) { \
  std::cout << "\nTesting " << N << std::endl; \
  T t=0; \
  if (!test_basic(t)) { \
    std::cout << "Error." << std::endl; return 1; \
  } \
}

int main()
{
  // builtin NTs
  TESTIT(int, "int")
  TESTIT(long int, "long int")
  TESTIT(short int, "short int")
  TESTIT(unsigned int, "unsigned int")
  TESTIT(unsigned long int, "unsigned long int")
  TESTIT(unsigned short int, "unsigned short int")
#ifdef LONG_LONG
  TESTIT(long long, "long long")
  TESTIT(unsigned long long, "unsigned long long")
#endif // LONG_LONG
  TESTIT(float, "float")
  TESTIT(double, "double")

  // CGAL number types
  TESTIT(CGAL::Fixed_precision_nt, "Fixed_precision_nt")
  TESTIT(CGAL::MP_Float, "MP_Float")
  TESTIT(CGAL::Quotient<int>, "Quotient<int>")
  TESTIT(QMPF, "Quotient<MP_Float>")
  TESTIT(CGAL::Lazy_exact_nt<QMPF>, "Lazy_exact_nt<Quotient<MP_Float> >")
#ifndef CGAL_CFG_MATCHING_BUG_2
#ifdef CGAL_USE_GMP
  TESTIT(FEIG, "Filtered_exact<int, Gmpz>");
  TESTIT(FEDG, "Filtered_exact<double, Gmpz>");
#endif
  TESTIT(FEDQ, "Filtered_exact<double, Quotient<MP_Float> >");
#endif // CGAL_CFG_MATCHING_BUG_2

  // GMP based NTs
#ifdef CGAL_USE_GMP
  TESTIT(CGAL::Gmpz, "Gmpz")
#endif // CGAL_USE_GMP

  // LEDA based NTs
#ifdef CGAL_USE_LEDA
  TESTIT(leda_integer, "leda_integer")
  TESTIT(leda_rational, "leda_rational")
  TESTIT(leda_bigfloat, "leda_bigfloat")
  TESTIT(leda_real, "leda_real")
#ifndef CGAL_CFG_MATCHING_BUG_2
  TESTIT(FEDR, "Filtered_exact<double, leda_real>");
#endif // CGAL_CFG_MATCHING_BUG_2
#endif // CGAL_USE_LEDA


  return 0;
}
