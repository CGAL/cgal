#include <CGAL/basic.h>
#include <CGAL/_test_utilities.h>

#include <CGAL/Quotient.h> 
#include <CGAL/MP_Float.h> 
#include <CGAL/Lazy_exact_nt.h> 
#include <CGAL/Fixed_precision_nt.h> 

#ifndef CGAL_CFG_MATCHING_BUG_2
#include <CGAL/Filtered_exact.h> 
#endif // CGAL_CFG_MATCHING_BUG_2

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#endif // CGAL_USE_GMP

#ifdef CGAL_USE_GMPXX
#include <CGAL/gmpxx.h>
#endif

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_Expr.h>
#endif

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

#define TESTIT(T,N) { \
  std::cout << "\nTesting " << N << std::endl; \
  T t=0; \
  if (!CGAL::test_utilities(t)) { \
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
  TESTIT(CGAL::Gmpq, "Gmpq")
#endif // CGAL_USE_GMP
#ifdef CGAL_USE_GMPXX
  TESTIT(mpz_class, "mpz_class")
  TESTIT(mpq_class, "mpq_class")
  TESTIT(mpf_class, "mpf_class")
#endif

  // CORE
#ifdef CGAL_USE_CORE
  TESTIT(CORE::Expr, "CORE::Expr")
#endif

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
