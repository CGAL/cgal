#include <CGAL/basic.h>
#include <CGAL/kernel_basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/predicates_on_ftC2.h>
// #include <CGAL/Gmpz.h>
#include <CGAL/leda_real.h>
// #include <CGAL/leda_bigfloat.h>
#include <CGAL/Filter.h>
// #define CGAL_IA_NO_WARNINGS
// #include <CGAL/Interval_arithmetic.h>

// Don't be stupid, CGAL_Gmpz can only store integers !!!
// typedef CGAL_Filtering<double, CGAL_Gmpz> NT;
typedef CGAL_Filtering<double, leda_real> NT;
// typedef CGAL_Filtering<double, leda_bigfloat> NT;
// typedef double NT;

// The predicate code is 100% generic.
// The #inclusion should be done automatically.
// #include <CGAL/Filter/predicates_on_ftC2.h>

// The program code is 100% generic/templated, as usual.
int main()
{
  NT px, py, la, lb, lc;
  px=1; py=2; la=3; lb=4; lc=5;
  cout << "Result 1st test: " << (int)CGAL_compare_y_at_xC2(px, py, la, lb, lc);
  cout << " ( == 1 )\n";
  px=1.1; py=1.7; la=1.3; lb=1.5; lc=-3.98;
  cout << "Result 2nd test: " << (int)CGAL_compare_y_at_xC2(px, py, la, lb, lc);
  cout << " ( == 0 ) (not sure, it depends of the first approx)\n";
}
