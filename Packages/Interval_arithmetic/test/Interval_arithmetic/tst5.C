#include <CGAL/basic.h>
#include <CGAL/kernel_basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/predicates_on_ftC2.h>
#include <CGAL/Gmpz.h>
#include <CGAL/leda_real.h>
#include <CGAL/Filter.h>
#define CGAL_IA_NO_WARNINGS
#include <CGAL/Interval_arithmetic.h>

typedef CGAL_Filtering<double, CGAL_Gmpz> NT;
// typedef CGAL_Filtering<double, leda_real> NT;
// typedef double NT;

// The predicate code is 100% generic.
#include <CGAL/filtered_predicates_on_ftC2.h>

// The program code is 100% generic/templated.
int main()
{
  NT px, py, la, lb, lc;
  px=1; py=2; la=3; lb=4; lc=5;
  cout << "Result 1st test: " << CGAL_compare_y_at_xC2(px, py, la, lb, lc);
  cout << endl;
  px=1.1; py=1.7; la=1.3; lb=1.5; lc=-3.98;
  cout << "Result 2nd test: " << (int)CGAL_compare_y_at_xC2(px, py, la, lb, lc);
  cout << endl;
}
