
// Test program for the Static Adaptatif filter.

#define CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

#include <CGAL/basic.h>
#include <iostream>

// Workaround for crappy compilers.
#ifdef CGAL_CFG_MATCHING_BUG_2
#define CGAL_IA_CT double
#define CGAL_IA_ET CGAL::MP_Float
#define CGAL_IA_CACHE No_Filter_Cache
#endif

#include <CGAL/MP_Float.h>
#include <CGAL/Filtered_exact.h>

#include <CGAL/predicates/kernel_ftC2.h>

typedef CGAL::Filtered_exact<double, CGAL::MP_Float, CGAL::Static> NT;
// typedef CGAL::Filtered_exact<double, CGAL::Filtered_exact<double, double>, CGAL::Static> NT;

int main()
{
  using std::cout;
  using std::endl;

  double px=0; double py=0;
  double qx=1; double qy=0;
  double rx=0; double ry=1;

  cout << "Call of orientationC2(";
  cout << px<<","<<py<<","<<qx<<","<<qy<<","<<rx<<","<<ry<<")\n";
  cout << "double : " << (int)CGAL::orientationC2(px,py,qx,qy,rx,ry);
  cout << endl;
  cout << "(SAF)  : ";
  cout << (int)CGAL::orientationC2(NT(px), NT(py), NT(qx), NT(qy), NT(rx), NT(ry));
  cout << endl;
  return 0;
}
