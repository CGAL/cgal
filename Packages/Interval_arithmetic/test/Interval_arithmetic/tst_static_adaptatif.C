
// Test program for the Static Adaptatif filter.

#define CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

#include <CGAL/basic.h>
#include <iostream>

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

// The number types.
#include <CGAL/Filtered_exact.h>
#include <CGAL/double.h>
// #include <CGAL/leda_real.h>

// The template predicate.
#include <CGAL/predicates/kernel_ftC2.h>

// typedef CGAL::Filtered_exact<double, leda_real, CGAL::Static> NT;
typedef CGAL::Filtered_exact<double, CGAL::Filtered_exact<double, double>, CGAL::Static> NT;

int main()
{
  using std::cout;
  using std::endl;
  using CGAL::orientationC2;

  double px=0; double py=0;
  double qx=1; double qy=0;
  double rx=0; double ry=1;

  cout << "Call of orientationC2(";
  cout << px<<","<<py<<","<<qx<<","<<qy<<","<<rx<<","<<ry<<")\n";
  cout << "double : " << (int)orientationC2(px,py,qx,qy,rx,ry);
  cout << endl;
  cout << "(SAF)  : ";
  cout << (int)orientationC2(NT(px), NT(py), NT(qx), NT(qy), NT(rx), NT(ry));
  cout << endl;
  return 0;
}
