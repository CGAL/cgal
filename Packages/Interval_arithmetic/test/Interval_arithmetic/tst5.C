#include <CGAL/config.h>
#include <CGAL/basic.h>

// Workaround for MipsPro.
#ifdef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
#define CGAL_IA_CT double
#define CGAL_IA_ET leda_real
#define CGAL_IA_CACHE No_Filter_Cache
#endif

#include <CGAL/misc.h>
// #include <CGAL/kernel_basic.h>
// #include <CGAL/number_utils.h>
#include <CGAL/Timer.h>

#include <CGAL/Quotient.h>
#include <CGAL/Interval_arithmetic.h>

#include <CGAL/predicates_on_ftC2.h>
// #include <CGAL/predicates_on_ftC3.h>

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

// #define CGAL_IA_NO_WARNINGS
// #include <CGAL/Interval_arithmetic.h>
// #include <CGAL/Fixed_precision_nt.h>
#include <CGAL/Double_eps.h>

#define CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
#include <CGAL/Arithmetic_filter.h>

using namespace CGAL;

// Please pay attention to the workaround for MipsPro at the top of the file.

// Don't be stupid, Gmpz can only store integers !!!
// typedef Filtered_exact<double, Gmpz> NT;
// typedef Filtered_exact<leda_real, leda_real> NT;
// typedef Filtered_exact<double, leda_bigfloat> NT;
#ifdef CGAL_USE_LEDA
typedef Filtered_exact<double, leda_real> NT;
#else
typedef Filtered_exact<double, Gmpz> NT; // Should be rationnals
#endif // CGAL_USE_LEDA
// typedef Filtered_exact<double, leda_real, Filter_Cache> NT;
// typedef Filtered_exact<double, leda_rational> NT;
// typedef Filtered_exact<float, leda_real> NT;
// typedef Filtered_exact<unsigned short int, leda_real> NT;
// typedef Filtered_exact<int, leda_real> NT;
// typedef Fixed_precision_nt NT;
// typedef double NT;
// typedef Gmpz NT;
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
  set_eps(0.0001);
  bench();
  return test();
}

void bench()
{
  const int loops = 1000000;
  int i;
  Timer t;
  double dt;
  Comparison_result result;
  NT px, py, la, lb, lc;

  px=1; py=2.0/3; la=3.0/5; lb=4.0/7; lc=5.0/6;
  px = py;
  px = NT(py);
  // NT ppx(py);
  dt = t.time(); t.start();
  for (i=0; i<loops; i++)
    result = compare_y_at_xC2(px, py, la, lb, lc);
  t.stop();
  std::cout << (int) result << "\t" << t.time()-dt << std::endl;
}

// The program code is 100% generic/templated, as usual.
int test()
{
  NT px, py, la, lb, lc;
  NT a (1);
#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
  Filtered_exact< Quotient<Gmpz>, Quotient<Gmpz> > qq (3,5);
  std::cout << (int) compare_y_at_xC2(qq,qq,qq,qq,qq);
#ifdef CGAL_USE_LEDA
  Filtered_exact< Quotient<int>, leda_rational> ii (3,2);
  Filtered_exact< Quotient<leda_integer>, Quotient<leda_integer> > jj (4,5);
  Interval_nt nt = jj.interval();
  std::cout << nt << std::endl;
#endif // CGAL_USE_LEDA
#endif // CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
  px=1; py=2; la=3; lb=4; lc=5;
  std::cout << "Result 1st test: " << (int)compare_y_at_xC2(px, py, la, lb, lc);
  std::cout << " ( == 1 )\n";
  px=1.1; py=1.7; la=1.3; lb=1.5; lc=-3.98;
  std::cout << "Result 2nd test: " << (int)compare_y_at_xC2(px, py, la, lb, lc);
  std::cout << " ( == 0 ) (not sure, it depends of the first approx)\n";
  std::cout << a << std::endl;
  return 0;
}
