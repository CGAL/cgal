// Test program for Lazy_exact_nt<>.

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Arithmetic_filter.h>
#include <CGAL/Lazy_exact_nt.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_real.h>
#else
typedef double leda_real; // :)
#endif

#ifdef __GNUG__
#define UNUSED __attribute__((unused))
#else
#define UNUSED
#endif

/*
  Compilation with GCC-2.95:
  g++ -W -Wall -Winline -fno-exceptions -fomit-frame-pointer -O3 -S Lazy_exact_nt.C
  the .s gives 3 __builtin_delete with GCC-2.95.
  GCC-2.96 does NULL propagation, and it works ;-)
*/

// typedef CGAL::Lazy_exact_nt<int> NT;
// typedef CGAL::Lazy_exact_nt<leda_real> NT;
typedef CGAL::Filtered_exact<CGAL::Lazy_exact_nt<leda_real>, leda_real> NT;
// typedef leda_real NT;

typedef CGAL::Cartesian<NT> Kernel;
typedef Kernel::Point_2     Point;

using std::cout;
using std::endl;

void predicats()
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
  // NT UNUSED a;
  // NT UNUSED b(a);
  // NT UNUSED c = b;
  NT UNUSED d (1.0);
  NT UNUSED e = d + d;
  NT UNUSED z = min(e,d);
  cout << e/NT(3) << endl;
  // NT UNUSED f = abs(NT(0));
  cout << "sign(3*(1/3)-1) = " << CGAL_NTS sign(NT(3)*(NT(1)/NT(3))-NT(1)) << endl;
  cout << "sign(sqrt(2)^2-2) = " << CGAL_NTS sign(CGAL_NTS square(CGAL_NTS sqrt(NT(2)))-NT(2)) << endl;
  // cout << "sign(sqrt(2)) = " << CGAL_NTS sign(CGAL_NTS sqrt(NT(2))) << endl;
  // cout << "sign(square(2)) = " << CGAL_NTS sign(CGAL_NTS square(NT(2))) << endl;
  // bool UNUSED pipo = e < d;
  cout << "sizeof(rep) = " << sizeof(CGAL::Lazy_exact_rep<int>) << endl;
  cout << "sizeof(cst) = " << sizeof(CGAL::Lazy_exact_Cst<int>) << endl;
  cout << "sizeof(abs) = " << sizeof(CGAL::Lazy_exact_Abs<int>) << endl;
  cout << "sizeof(add) = " << sizeof(CGAL::Lazy_exact_Add<int>) << endl;
  cout << "sizeof(leda_real) = " << sizeof(leda_real) << endl;
  cout << "sizeof(Lazy_exact_nt) = " << sizeof(CGAL::Lazy_exact_nt<int>) << endl;
  predicats();
  return 0;
}

