
#include <CGAL/basic.h>
#include <CGAL/config.h>
// #include <CGAL/kernel_basic.h>
// #include <CGAL/number_utils.h>
#include <CGAL/Timer.h>

#include <CGAL/predicates_on_ftC2.h>
// #include <CGAL/Arithmetic_filter/predicates_on_ftC2.h>

#define CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
#include <CGAL/Arithmetic_filter.h>

#include <CGAL/Gmpz.h>
#include <CGAL/leda_real.h>
#include <CGAL/double.h>
// #include <CGAL/leda_bigfloat.h>

// #define CGAL_IA_NO_WARNINGS
// #include <CGAL/Interval_arithmetic.h>
// #include <CGAL/Fixed.h>
#include <CGAL/Double_eps.h>
#include <CGAL/predicates_on_ftC3.h>
#include <CGAL/Arithmetic_filter/predicates_on_ftC3.h>
#include <CGAL/predicates_on_rtH2.h>
#include <CGAL/Arithmetic_filter/predicates_on_rtH2.h>

// Don't be stupid, CGAL_Gmpz can only store integers !!!
// typedef CGAL_Filtered_exact<double, CGAL_Gmpz> NT;
// typedef CGAL_Filtered_exact<leda_real, leda_real> NT;
// typedef CGAL_Filtered_exact<double, leda_bigfloat> NT;
typedef CGAL_Filtered_exact<double, leda_real> NT;
// typedef CGAL_Filtered_exact<float, leda_real> NT;
// typedef CGAL_Filtered_exact<unsigned short int, leda_real> NT;
// typedef CGAL_Filtered_exact<int, leda_real> NT;
// typedef Fixed NT;
// typedef CGAL_Gmpz NT;
// typedef double NT;
// typedef leda_real NT;
// typedef leda_rational NT;
// typedef CGAL_Double_eps NT;
// typedef CGAL_Interval_nt NT;
// typedef CGAL_Interval_nt_advanced NT;

int test(void);
void bench(void);

int main()
{
  // CGAL_Interval_nt ia (2);
  // double d = CGAL_to_double(ia);
  // cin >> ia;
  // cout << ia << endl;
  CGAL_set_eps(0.0001);
  bench();
  return test();
}

void bench()
{
  const int loops = 1000000;
  int i;
  CGAL_Timer t;
  double dt;
  CGAL_Comparison_result result;
  NT px, py, la, lb, lc;

  px=1; py=2; la=3; lb=4; lc=5;
  px = py;
  px = NT(py);
  NT ppx(py);
  dt = t.time(); t.start();
  for (i=0; i<loops; i++)
    result = CGAL_compare_y_at_xC2(px, py, la, lb, lc);
  t.stop();
  cout << (int) result << "\t" << t.time()-dt << endl;
}

// The program code is 100% generic/templated, as usual.
int test()
{
  NT px, py, la, lb, lc;
  px=1; py=2; la=3; lb=4; lc=5;
  cout << "Result 1st test: " << (int)CGAL_compare_y_at_xC2(px, py, la, lb, lc);
  cout << " ( == 1 )\n";
  px=1.1; py=1.7; la=1.3; lb=1.5; lc=-3.98;
  cout << "Result 2nd test: " << (int)CGAL_compare_y_at_xC2(px, py, la, lb, lc);
  cout << " ( == 0 ) (not sure, it depends of the first approx)\n";
  return 0;
}
