#include <CGAL/basic.h>

// Workaround for crappy compilers.
#ifdef CGAL_CFG_MATCHING_BUG_2
#define CGAL_IA_CT double
#ifdef CGAL_USE_LEDA
#define CGAL_IA_ET leda_real
#else
#define CGAL_IA_ET double
#endif
#define CGAL_IA_CACHE No_Filter_Cache
#endif

#include <CGAL/misc.h>
#include <CGAL/number_utils.h>
#include <CGAL/Timer.h>

#include <CGAL/Quotient.h>
#include <CGAL/Interval_arithmetic.h>

#include <CGAL/predicates/kernel_ftC2.h>
// #include <CGAL/predicates/kernel_ftC3.h>

#include <CGAL/double.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif

#ifdef CGAL_USE_LEDA
// #include <CGAL/leda_bigfloat.h>
// #include <CGAL/leda_rational.h>
#include <CGAL/leda_real.h>
#include <CGAL/leda_integer.h>
#endif // CGAL_USE_LEDA

// #include <CGAL/Interval_arithmetic.h>
// #include <CGAL/Fixed_precision_nt.h>
#include <CGAL/Double_eps.h>

#define CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
#include <CGAL/Arithmetic_filter.h>

// Please pay attention to the workaround at the top of the file.

// Don't be stupid, Gmpz can only store integers !!!
// typedef CGAL::Filtered_exact<double, CGAL::Gmpz> NT;
// typedef CGAL::Filtered_exact<leda_real, leda_real> NT;
// typedef CGAL::Filtered_exact<double, leda_bigfloat> NT;
#ifdef CGAL_USE_LEDA
typedef CGAL::Filtered_exact<double, leda_real> NT;
#elif defined(CGAL_USE_GMP)
typedef CGAL::Filtered_exact<double, CGAL::Gmpz> NT; // Should be exact rationnals
#else
typedef CGAL::Filtered_exact<double, double> NT; // Should be exact rationnals
#endif // CGAL_USE_LEDA
// typedef CGAL::Filtered_exact<double, leda_real, Filter_Cache> NT;
// typedef CGAL::Filtered_exact<double, leda_rational> NT;
// typedef CGAL::Filtered_exact<float, leda_real> NT;
// typedef CGAL::Filtered_exact<unsigned short int, leda_real> NT;
// typedef CGAL::Filtered_exact<int, leda_real> NT;
// typedef Fixed_precision_nt NT;
// typedef double NT;
// typedef CGAL::Gmpz NT;
// typedef leda_real NT;
// typedef leda_rational NT;
// typedef Double_eps NT;
// typedef Interval_nt NT;
// typedef Interval_nt_advanced NT;

int test(void);
void bench(void);

// Now try the higher level predicates (C2/C3/H2/H3).
// Well, that is make a test-suite.
int main()
{
  // Interval_nt ia (2);
  // double d = to_double(ia);
  // cin >> ia;
  // std::cout << ia << std::endl;
  CGAL::set_eps(0.0001);
  bench();
  return test();
}

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

// The program code is 100% generic/templated, as usual.
int test()
{
  NT px, py, la, lb, lc;
  NT a (1);
  // a = CGAL::abs(a);
#ifndef CGAL_CFG_MATCHING_BUG_2
#ifdef CGAL_USE_GMP
  CGAL::Filtered_exact< CGAL::Quotient<CGAL::Gmpz>, CGAL::Quotient<CGAL::Gmpz> > qq (3,5);
  std::cout << (int) CGAL::compare_y_at_xC2(qq,qq,qq,qq,qq);
#endif
#ifdef CGAL_USE_LEDA
  CGAL::Filtered_exact< CGAL::Quotient<int>, leda_rational> ii (3,2);
  CGAL::Filtered_exact< CGAL::Quotient<leda_integer>, CGAL::Quotient<leda_integer> > jj (4,5);

  CGAL::FPU_CW_t backup = CGAL::FPU_get_and_set_cw(CGAL::FPU_cw_up);
  CGAL::Interval_nt nt = jj.interval();
  CGAL::FPU_set_cw(backup);

  std::cout << nt << std::endl;
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
