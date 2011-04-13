/*
 * This test file is not extensive enough.
 */

#include <CGAL/basic.h>

// Workaround for buggy compilers.
#ifdef CGAL_CFG_MATCHING_BUG_2
#define CGAL_IA_CT double
#define CGAL_IA_PROTECTED true
#define CGAL_IA_CACHE No_Filter_Cache
#ifdef CGAL_USE_LEDA
#  define CGAL_IA_ET leda_real
#elif defined CGAL_USE_GMP
#  define CGAL_IA_ET CGAL::Quotient<CGAL::Gmpz>
#else
#  define CGAL_IA_ET CGAL::MP_Float
#endif
#endif

#include <CGAL/Timer.h>

#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_real.h>
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#endif

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif

#include <CGAL/predicates/kernel_ftC2.h>
// #include <CGAL/predicates/kernel_ftC3.h>

#define CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
#include <CGAL/Filtered_exact.h>

// PLEASE PAY ATTENTION TO THE WORKAROUND AT THE TOP OF THE FILE !!!
#ifdef CGAL_USE_LEDA
typedef CGAL::Filtered_exact<double, leda_real> NT;
// #elif defined CGAL_USE_GMP
// typedef CGAL::Filtered_exact<double, CGAL::Quotient<CGAL::Gmpz> > NT;
#else
typedef CGAL::Filtered_exact<double, CGAL::Quotient<CGAL::MP_Float> > NT;
#endif


void bench()
{
  const int loops = 1000000;
  int i;
  CGAL::Timer t;
  double dt;
  CGAL::Comparison_result result;
  NT px, py, la, lb, lc;

  px=1; py=2.0/3; la=3.0/5; lb=4.0/7; lc=5.0/6;
  px = py;
  px = NT(py);
  // NT ppx(py);
  dt = t.time(); t.start();
  for (i=0; i<loops; i++)
    result = CGAL::compare_y_at_xC2(px, py, la, lb, lc);
  t.stop();
  std::cout << (int) result << "\t" << t.time()-dt << std::endl;
}

#ifndef CGAL_CFG_MATCHING_BUG_2
namespace CGAL {
template <class NT>
NT
my_abs (const NT &n)
{
  return CGAL_NTS abs(n);
}
} // namespace CGAL
#endif

// The program code is 100% generic/templated, as usual.
int test()
{
  NT px, py, la, lb, lc;
  NT a (1);
#ifndef CGAL_CFG_MATCHING_BUG_2
  a = CGAL::my_abs(a);
#ifdef CGAL_USE_GMP
  CGAL::Filtered_exact< CGAL::Quotient<CGAL::Gmpz>, CGAL::Quotient<CGAL::Gmpz> > qq (3,5);
  std::cout << (int) CGAL::compare_y_at_xC2(qq,qq,qq,qq,qq);
#endif
#ifdef CGAL_USE_LEDA
  CGAL::Filtered_exact< CGAL::Quotient<int>, leda_rational> ii (3,2);
  CGAL::Filtered_exact< CGAL::Quotient<leda_integer>, CGAL::Quotient<leda_integer> > jj (4,5);

  CGAL::FPU_CW_t backup = CGAL::FPU_get_and_set_cw(CGAL_FE_UPWARD);
  CGAL::Interval_nt_advanced nt = jj.interval();
  CGAL::Interval_nt_advanced nt2 = ii.interval();
  CGAL::FPU_set_cw(backup);

  std::cout << nt << std::endl;
  std::cout << nt2 << std::endl;
#endif // CGAL_USE_LEDA
#endif // CGAL_CFG_MATCHING_BUG_2
  px=1; py=2; la=3; lb=4; lc=5;
  std::cout << "Result 1st test: " << (int)CGAL::compare_y_at_xC2(px, py, la, lb, lc);
  std::cout << " ( == 1 )\n";
  px=1.1; py=1.7; la=1.3; lb=1.5; lc=-3.98;
  std::cout << "Result 2nd test: " << (int)CGAL::compare_y_at_xC2(px, py, la, lb, lc);
  std::cout << " ( == 0 ) (not sure, it depends of the first approx)\n";
  std::cout << a << std::endl;
  return 0;
}

int main()
{
  bench();
  return test();
}
