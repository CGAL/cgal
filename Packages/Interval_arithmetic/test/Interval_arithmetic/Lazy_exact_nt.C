
// Test program for Lazy_exact_nt<>.

#define CGAL_NO_ASSERTIONS

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Lazy_exact_nt.h>

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

typedef CGAL::Lazy_exact_nt<int> NT;

int main ()
{
  NT UNUSED a;
  NT UNUSED b(a);
  NT UNUSED c = b;
  NT UNUSED d = 1.0;
  NT UNUSED e = d + d;
  bool UNUSED pipo = e < d;
  return 0;
}

