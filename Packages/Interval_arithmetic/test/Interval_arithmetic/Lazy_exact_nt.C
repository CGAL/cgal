
// Test program for Lazy_exact_nt<>.

#define CGAL_NO_ASSERTIONS

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/leda_real.h>

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
typedef CGAL::Lazy_exact_nt<leda_real> NT;
// typedef leda_real NT;

int main ()
{
  using std::cout;
  using std::endl;
  cout.precision(20);
  // NT UNUSED a;
  // NT UNUSED b(a);
  // NT UNUSED c = b;
  NT UNUSED d = 1.0;
  NT UNUSED e = d + d;
  NT UNUSED z = min(e,d);
  cout << e/3 << endl;
  // NT UNUSED f = abs(NT(0));
  cout << "sign(3*(1/3)-1) = " << CGAL_NTS::sign(NT(3)*(NT(1)/NT(3))-NT(1)) << endl;
  cout << "sign(sqrt(2)^2-2) = " << CGAL_NTS::sign(CGAL_NTS::square(CGAL_NTS::sqrt(NT(2)))-NT(2)) << endl;
  // cout << "sign(sqrt(2)) = " << CGAL_NTS::sign(CGAL_NTS::sqrt(NT(2))) << endl;
  // cout << "sign(square(2)) = " << CGAL_NTS::sign(CGAL_NTS::square(NT(2))) << endl;
  // bool UNUSED pipo = e < d;
  cout << "sizeof(rep) = " << sizeof(CGAL::Lazy_exact_nt_rep<int>) << endl;
  cout << "sizeof(cst) = " << sizeof(CGAL::Lazy_exact_nt_Cst<int>) << endl;
  cout << "sizeof(abs) = " << sizeof(CGAL::Lazy_exact_nt_Abs<int>) << endl;
  cout << "sizeof(add) = " << sizeof(CGAL::Lazy_exact_nt_Add<int>) << endl;
  cout << "sizeof(leda_real) = " << sizeof(leda_real) << endl;
  cout << "sizeof(Lazy_exact_nt) = " << sizeof(CGAL::Lazy_exact_nt<int>) << endl;
  return 0;
}

