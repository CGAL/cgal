#include <CGAL/basic.h>
#include <CGAL/Quotient.h> 
#include <CGAL/MP_Float.h> 
#include <CGAL/Lazy_exact_nt.h> 
#include <CGAL/Fixed_precision_nt.h> 
#include <CGAL/Filtered_exact.h> 

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif // CGAL_USE_GMP

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>
#endif // CGAL_USE_LEDA

template < class NT >
bool
test_sqrt(NT, CGAL::Tag_true)
{
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
  NT eleven(11);
  NT three(3);
  std::cerr << "  div()" << std::endl;
  if (CGAL_NTS div(eleven, three) != three) return false;
  NT x(2737);
  NT y(5083);
  std::cerr << "  gcd()" << std::endl;
  if (CGAL_NTS gcd(x, y) != NT(391)) return false;

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
  std::cerr << "  is_zero()" << std::endl;
  if (!CGAL_NTS is_zero(zero)) return false;
  std::cerr << "  is_one()" << std::endl;
  if (!CGAL_NTS is_one(one)) return false;
  std::cerr << "  is_positive()" << std::endl;
  if (!CGAL_NTS is_positive(one)) return false;
  std::cerr << "  is_negative()" << std::endl;
  if (CGAL_NTS is_negative(one)) return false;
  std::cerr << "  sign()" << std::endl;
  if (CGAL_NTS sign(one) != CGAL::POSITIVE) return false;
  if (CGAL_NTS sign(zero) != CGAL::ZERO) return false;
  std::cerr << "  abs()" << std::endl;
  if (CGAL_NTS abs(-one) != one &&
      CGAL_NTS abs(-one) != -one) // unsigned types :-)
    return false;
  std::cerr << "  compare()" << std::endl;
  if (CGAL_NTS compare(one, zero) != CGAL::LARGER) return false;
  if (CGAL_NTS compare(one, one) != CGAL::EQUAL) return false;
  if (CGAL_NTS compare(zero, one) != CGAL::SMALLER) return false;
  std::cerr << "  square()" << std::endl;
  NT two(2);
  NT four = CGAL_NTS square(two);

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
  TESTIT(CGAL::MP_Float, "CGAL::MP_Float")
  TESTIT(CGAL::Quotient<int>, "CGAL::Quotient<int>")
    //  TESTIT(CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> >,
    //	 "CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> >")

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
#endif // CGAL_USE_LEDA


  return 0;
}
