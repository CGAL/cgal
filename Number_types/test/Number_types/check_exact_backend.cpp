#include <CGAL/config.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Exact_rational.h>

#include <type_traits>

using namespace CGAL::internal;

int main()
{
  static_assert(
#ifdef CGAL_USE_GMP
    ( Default_exact_nt_backend!=GMP_BACKEND || (std::is_same_v<CGAL::Exact_integer,CGAL::Gmpz> && std::is_same_v<CGAL::Exact_rational,CGAL::Gmpq>) ) &&
#endif
#ifdef CGAL_USE_GMPXX
    ( Default_exact_nt_backend!=GMPXX_BACKEND || (std::is_same_v<CGAL::Exact_integer,mpz_class> && std::is_same_v<CGAL::Exact_rational,mpq_class>) ) &&
#endif
#if defined(CGAL_USE_BOOST_MP) && defined(CGAL_USE_GMP)
    ( Default_exact_nt_backend!=BOOST_GMP_BACKEND || (std::is_same_v<CGAL::Exact_integer,boost::multiprecision::mpz_int> && std::is_same_v<CGAL::Exact_rational,boost::multiprecision::mpq_rational>) ) &&
#endif
#if defined(CGAL_USE_BOOST_MP)
    ( Default_exact_nt_backend!=BOOST_BACKEND || (std::is_same_v<CGAL::Exact_integer,boost::multiprecision::cpp_int> && std::is_same_v<CGAL::Exact_rational,boost::multiprecision::cpp_rational>) ) &&
#endif
#if defined(CGAL_USE_LEDA)
    ( Default_exact_nt_backend!=LEDA_BACKEND || (std::is_same_v<CGAL::Exact_integer,leda_integer> && std::is_same_v<CGAL::Exact_rational,leda_rational>) ) &&
#endif
    ( Default_exact_nt_backend!=MP_FLOAT_BACKEND || (std::is_same_v<CGAL::Exact_integer,CGAL::MP_Float> && std::is_same_v<CGAL::Exact_rational,CGAL::Quotient<CGAL::MP_Float>>) )
  );

  std::cout << "Exact backend is " << exact_nt_backend_string() << "\n";
#ifdef CGAL_USE_CORE
#ifdef CGAL_CORE_USE_GMP_BACKEND
  std::cout << "CGAL_CORE_USE_GMP_BACKEND is defined, using gmp backend in BigInt and BigRat\n";
#else
  std::cout << "CGAL_CORE_USE_GMP_BACKEND is NOT defined, using boost-mp backend in BigInt and BigRat\n";
#endif
#else
  std::cout << "CGAL_USE_CORE is not defined\n";
#endif
}
