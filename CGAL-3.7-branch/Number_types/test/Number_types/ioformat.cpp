#include <CGAL/basic.h>

#include <CGAL/Quotient.h> 
#include <CGAL/MP_Float.h> 
#include <CGAL/Lazy_exact_nt.h> 
#include <CGAL/Interval_nt.h> 

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/Gmpq.h>
#endif

#ifdef CGAL_USE_GMPXX
#include <CGAL/gmpxx.h>
#endif

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/CORE_BigFloat.h>
#include <CGAL/CORE_Expr.h>
#endif

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>
#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_LONG_LONG
#  include <CGAL/long_long.h>
#endif

#include <CGAL/Number_type_checker.h>
#include <cassert>
#include <CGAL/_test_utilities.h>

typedef CGAL::Quotient<CGAL::MP_Float>            QMPF;

#define TESTIT(NT,N)                                                    \
    {                                                                   \
        std::cout << "\nTesting ioformat: " << N << std::endl;          \
        NT tmp2(0), tmp1(13);                                           \
                                                                        \
        std::ostringstream os;                                          \
        os << ::CGAL::oformat(tmp1);                                    \
        std::istringstream is(os.str());                                \
        is >> ::CGAL::iformat(tmp2);                                    \
        assert( tmp1 == tmp2 );                                         \
    }
    

int main()
{
  
  // builtin NTs
  TESTIT(int, "int")
  TESTIT(long int, "long int")
  TESTIT(short int, "short int")
  // Unsigned types are not appropriate for many things...
  // TESTIT(unsigned int, "unsigned int")
  // TESTIT(unsigned long int, "unsigned long int")
  // TESTIT(unsigned short int, "unsigned short int")
#ifdef CGAL_USE_LONG_LONG
  TESTIT(long long, "long long")
  // TESTIT(unsigned long long, "unsigned long long")
#endif
  TESTIT(float, "float")
  TESTIT(double, "double")
  TESTIT(long double, "long double")

  // CGAL number types
  //TESTIT(CGAL::MP_Float, "MP_Float")
  TESTIT(CGAL::Quotient<int>, "Quotient<int>")
  TESTIT(QMPF, "Quotient<MP_Float>")
  TESTIT(CGAL::Lazy_exact_nt<QMPF>, "Lazy_exact_nt<Quotient<MP_Float> >")
  TESTIT(CGAL::Interval_nt<>, "Interval_nt<>")

  // GMP based NTs
#ifdef CGAL_USE_GMP
  TESTIT(CGAL::Gmpz, "Gmpz")
  TESTIT(CGAL::Gmpz, "Gmpzf")
  TESTIT(CGAL::MP_Float, "MP_Float")
  TESTIT(CGAL::Gmpq, "Gmpq")
#endif // CGAL_USE_GMP
#ifdef CGAL_USE_GMPXX
  TESTIT(mpz_class, "mpz_class")
  TESTIT(mpq_class, "mpq_class")
  // TESTIT(mpf_class, "mpf_class") // Not finished.
#endif

  // CORE
#ifdef CGAL_USE_CORE
      //bug in io for CORE. 
      TESTIT(CORE::BigInt, "CORE::BigInt")
      TESTIT(CORE::BigRat, "CORE::BigRat")
      TESTIT(CORE::BigFloat, "CORE::BigFloat")
      //TESTIT(CORE::Expr, "CORE::Expr")
#endif

      // LEDA based NTs
#ifdef CGAL_USE_LEDA
      TESTIT(leda_integer, "leda_integer")
      TESTIT(leda_rational, "leda_rational")
      TESTIT(leda_bigfloat, "leda_bigfloat")
      TESTIT(leda_real, "leda_real")
      typedef CGAL::Number_type_checker<leda_rational,leda_real> NT_checker;
      TESTIT(NT_checker, "NT_checker");
#endif // CGAL_USE_LEDA
      
  return 0;
}
