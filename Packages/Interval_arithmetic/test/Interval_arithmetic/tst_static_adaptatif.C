
// Test program for the Static Adaptatif filter.

#define CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

#include <CGAL/basic.h>
#include <math.h>
#include <iostream>

// The number types.
#include <CGAL/Arithmetic_filter.h>
#include <CGAL/double.h>
#include <CGAL/leda_real.h>
#include <CGAL/Static_adaptatif_filter.h>
#include <CGAL/Static_filter_error.h>
#include <CGAL/Restricted_double.h>

// The template predicate.
#include <CGAL/predicates_on_ftC2.h>
// The overloaded predicate.
#include <CGAL/Static_adaptatif_filter/predicates/sign_of_determinant.h>
#include <CGAL/Static_adaptatif_filter/predicates_on_ftC2.h>

int main()
{
  using std::cout;
  using std::endl;
  using CGAL::orientationC2;
  using CGAL::Static_adaptatif_filter;

  double px=0; double py=0;
  double qx=1; double qy=0;
  double rx=0; double ry=1;

  cout << "Call of orientationC2(";
  cout <<px<<","<<py<<","<<qx<<","<<qy<<","<<rx<<","<<ry<<")\n";
  cout << "double : " << (int)orientationC2(px,py,qx,qy,rx,ry) << endl;
  cout << "(SAF)  : " << (int)orientationC2(Static_adaptatif_filter(px),
  				Static_adaptatif_filter(py),
				Static_adaptatif_filter(qx),
				Static_adaptatif_filter(qy),
				Static_adaptatif_filter(rx),
				Static_adaptatif_filter(ry)) << endl;
  return 0;
}
