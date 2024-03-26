#include <CGAL/config.h>

#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Sqrt_extension.h>

#ifdef CGAL_USE_BOOST_MP
#include <CGAL/boost_mp.h>
#endif

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/Mpzf.h>
#include <CGAL/Gmpq.h>
#endif

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/CORE_BigFloat.h>
#include <CGAL/CORE_Expr.h>
#endif

#ifdef CGAL_USE_GMPXX
#include <CGAL/gmpxx.h>
#endif

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>
#endif // CGAL_USE_LEDA

#include <CGAL/long_long.h>

#include <CGAL/Number_type_checker.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/_test_utilities.h>

typedef CGAL::Quotient<CGAL::MP_Float>            QMPF;

#define TESTIT(T,N)                                             \
{                                                               \
    std::cout << "\nTesting " << N << std::endl;                \
    T t=0;                                                      \
    if (!CGAL::test_utilities(t)) {                             \
        std::cout << "Error." << std::endl; return 1;           \
    }                                                           \
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
  TESTIT(long long, "long long")
  // TESTIT(unsigned long long, "unsigned long long")
  TESTIT(float, "float")
  TESTIT(double, "double")
  TESTIT(long double, "long double")

  // CGAL number types
  //TESTIT(CGAL::MP_Float, "MP_Float") // CGAL::div(MP_Float, MP_Float) does not implement _integer_ division
  TESTIT(CGAL::Quotient<int>, "Quotient<int>")
  TESTIT(QMPF, "Quotient<MP_Float>")
  TESTIT(CGAL::Lazy_exact_nt<QMPF>, "Lazy_exact_nt<Quotient<MP_Float> >")
  TESTIT(CGAL::Interval_nt<>, "Interval_nt<>")

  // Boost.Multiprecision
#ifdef CGAL_USE_BOOST_MP
  TESTIT(boost::multiprecision::cpp_int, "cpp_int")
  TESTIT(boost::multiprecision::cpp_rational, "cpp_rational")
#endif

  // GMP based NTs
#ifdef CGAL_USE_GMP
  TESTIT(CGAL::Gmpz, "Gmpz")
  TESTIT(CGAL::Gmpzf, "Gmpzf")
# ifdef CGAL_HAS_MPZF
  TESTIT(CGAL::Mpzf, "Mpzf")
# endif
  TESTIT(CGAL::Gmpq, "Gmpq")
# ifdef CGAL_USE_BOOST_MP
  TESTIT(boost::multiprecision::mpz_int, "mpz_int")
  TESTIT(boost::multiprecision::mpq_rational, "mpq_rational")
# endif
#endif // CGAL_USE_GMP
#ifdef CGAL_USE_GMPXX
  TESTIT(mpz_class, "mpz_class")
  TESTIT(mpq_class, "mpq_class")
  // TESTIT(mpf_class, "mpf_class") // Not finished.
#endif

  // CORE
#ifdef CGAL_USE_CORE
      TESTIT(CORE::BigInt, "CORE::BigInt")
      TESTIT(CORE::BigRat, "CORE::BigRat")
      TESTIT(CORE::BigFloat, "CORE::BigFloat")
      TESTIT(CORE::Expr, "CORE::Expr")
      typedef CGAL::Number_type_checker<CORE::BigRat,CORE::Expr> NT_checker;
      TESTIT(NT_checker, "NT_checker");
#endif

      // LEDA based NTs
#ifdef CGAL_USE_LEDA
      TESTIT(leda_integer, "leda_integer")
      TESTIT(leda_rational, "leda_rational")
      TESTIT(leda_bigfloat, "leda_bigfloat")
      TESTIT(leda_real, "leda_real")
#endif // CGAL_USE_LEDA

       // TEST Sqrt_extension
#ifdef CGAL_USE_GMP
      typedef CGAL::Sqrt_extension<int,int> Ext_int;
      TESTIT(Ext_int     , "CGAL::Sqrt_extension<int,int>");
      typedef CGAL::Sqrt_extension<CGAL::Gmpz,CGAL::Gmpz> Ext_int_int;
      TESTIT(Ext_int_int , "CGAL::Sqrt_extension<CGAL::Gmpz,CGAL::Gmpz>");
      typedef CGAL::Sqrt_extension<CGAL::Gmpq,CGAL::Gmpz> Ext_rat_int;
      TESTIT(Ext_rat_int , "CGAL::Sqrt_extension<CGAL::Gmpq,CGAL::Gmpz>");
      typedef CGAL::Sqrt_extension<CGAL::Gmpq,CGAL::Gmpq> Ext_rat_rat;
      TESTIT(Ext_rat_rat , "CGAL::Sqrt_extension<CGAL::Gmpq,CGAL::Gmpq>");
#endif // CGAL_USE_GMP


  return 0;
}
