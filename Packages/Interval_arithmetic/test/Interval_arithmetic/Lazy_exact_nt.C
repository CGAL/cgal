// Test program for Lazy_exact_nt<>.

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Arithmetic_filter.h>
#include <CGAL/Lazy_exact_nt.h>

#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
typedef leda_real        Exact_NT;
#else
#  include <CGAL/MP_Float.h>
typedef CGAL::MP_Float   Exact_NT; // doesn't do exact sqrt() though
#endif

// typedef CGAL::Lazy_exact_nt<int> NT;
// typedef CGAL::Lazy_exact_nt<Exact_NT> NT;
typedef CGAL::Filtered_exact<CGAL::Lazy_exact_nt<Exact_NT>, Exact_NT> NT;
// typedef Exact_NT NT;

typedef CGAL::Cartesian<NT> Kernel;
typedef Kernel::Point_2     Point;

using std::cout;
using std::endl;

void predicates()
{
  Point A(NT(1.0)/NT(3),NT(2.0)/NT(3));
  Point B(NT(2.0)/NT(3),NT(3.0)/NT(3));
  Point C(NT(3.0)/NT(3),NT(4.0)/NT(3));
  Point D(NT(4.0)/NT(3),NT(3.0)/NT(3));
  cout << "A : " << A << endl;
  cout << "B : " << B << endl;
  cout << "C : " << C << endl;
  cout << (int) CGAL::orientation(A,B,C) << endl;
}

int main ()
{
  cout.precision(20);
  // NT a;
  // NT b(a);
  // NT c = b;
  NT d (1.0);
  NT e = d + d;
  NT z = min(e,d);
  (void) d; (void) e; (void) z; // Shut up warnings.
  cout << e/NT(3) << endl;
  // NT f = abs(NT(0));
  cout << "sign(3*(1/3)-1) = " << CGAL_NTS sign(NT(3)*(NT(1)/NT(3))-NT(1)) << endl;
  cout << "sign(sqrt(2)^2-2) = " << CGAL_NTS sign(CGAL_NTS square(CGAL_NTS sqrt(NT(2)))-NT(2)) << endl;
  // cout << "sign(sqrt(2)) = " << CGAL_NTS sign(CGAL_NTS sqrt(NT(2))) << endl;
  // cout << "sign(square(2)) = " << CGAL_NTS sign(CGAL_NTS square(NT(2))) << endl;
  // bool pipo = e < d;
  cout << "sizeof(rep) = " << sizeof(CGAL::Lazy_exact_rep<int>) << endl;
  cout << "sizeof(cst) = " << sizeof(CGAL::Lazy_exact_Cst<int>) << endl;
  cout << "sizeof(abs) = " << sizeof(CGAL::Lazy_exact_Abs<int>) << endl;
  cout << "sizeof(add) = " << sizeof(CGAL::Lazy_exact_Add<int>) << endl;
  cout << "sizeof(Exat_NT) = " << sizeof(Exact_NT) << endl;
  cout << "sizeof(Lazy_exact_nt) = " << sizeof(CGAL::Lazy_exact_nt<int>) << endl;
  predicates();
  return 0;
}

