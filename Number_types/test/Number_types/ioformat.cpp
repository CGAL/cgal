#include <CGAL/config.h>

#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Interval_nt.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/Mpzf.h>
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

#include <CGAL/long_long.h>

#include <CGAL/Number_type_checker.h>
#include <cassert>

#include <CGAL/disable_warnings.h>

#include <CGAL/_test_utilities.h>

typedef CGAL::Quotient<CGAL::MP_Float>            QMPF;

template <typename NT>
void test_it(const char* N, int value)
{
  std::cout << "\nTesting ioformat: " << N
            << " with " << value << std::endl;
  NT tmp2(0), tmp1 = static_cast<NT>(value);

  std::ostringstream os;
  os << ::CGAL::IO::oformat(tmp1);
  std::cout << os.str() << std::endl;
  std::istringstream is(os.str());
  is >> ::CGAL::IO::iformat(tmp2);
  assert( tmp1 == tmp2 );
}

template <typename NT>
void test_it(const char* N)
{
  test_it<NT>(N, 13);
  test_it<NT>(N, -27);
  test_it<NT>(N, 0);
}

int main()
{

  // builtin NTs
  test_it<int>("int");
  test_it<long int>("long int");
  test_it<short int>("short int");
  // Unsigned types are not appropriate for many things...
  // test_it<unsigned int>("unsigned int");
  // test_it<unsigned long int>("unsigned long int");
  // test_it<unsigned short int>("unsigned short int");
  test_it<long long>("long long");
  // test_it<unsigned long long>("unsigned long long");
  test_it<float>("float");
  test_it<double>("double");
  test_it<long double>("long double");

  // CGAL number types
  //test_it<CGAL::MP_Float>("MP_Float");
  test_it<CGAL::Quotient<int> >("Quotient<int>");
  test_it<QMPF>("Quotient<MP_Float>");
  test_it<CGAL::Lazy_exact_nt<QMPF> >("Lazy_exact_nt<Quotient<MP_Float> >");
  test_it<CGAL::Interval_nt<> >("Interval_nt<>");

  // GMP based NTs
#ifdef CGAL_USE_GMP
  test_it<CGAL::Gmpz>("Gmpz");
  test_it<CGAL::Gmpzf>("Gmpzf");
# ifdef CGAL_HAS_MPZF
  test_it<CGAL::Mpzf>("Mpzf");
# endif
  test_it<CGAL::MP_Float>("MP_Float");
  test_it<CGAL::Gmpq>("Gmpq");
#endif // CGAL_USE_GMP
#ifdef CGAL_USE_GMPXX
  test_it<mpz_class>("mpz_class");
  test_it<mpq_class>("mpq_class");
  // test_it<mpf_class>("mpf_class"); // Not finished.
#endif

  // CORE
#ifdef CGAL_USE_CORE
  //bug in io for CORE.
  test_it<CORE::BigInt>("CORE::BigInt");
  test_it<CORE::BigRat>("CORE::BigRat");
  test_it<CORE::BigFloat>("CORE::BigFloat");
  //test_it<CORE::Expr>("CORE::Expr");
#endif

  // LEDA based NTs
#ifdef CGAL_USE_LEDA
  static_assert(CGAL::Output_rep<leda_rational>::is_specialized == true);
  test_it<leda_integer>("leda_integer");
  test_it<leda_rational>("leda_rational");
  test_it<leda_bigfloat>("leda_bigfloat");
  test_it<leda_real>("leda_real");
  typedef CGAL::Number_type_checker<leda_rational,leda_real> NT_checker;
  test_it<NT_checker>("NT_checker");;
#endif // CGAL_USE_LEDA

  return 0;
}
