// Generic bench file for the IA package. include/bench_generic.C
//  $Revision$
//  $Date$
// Written by Sylvain Pion, 1997/1998.

// This file is included from tst[34].C, that do just a #define:
// #define TESTED_TYPE CGAL_Interval_nt_advanced // For tst3.C
// #define TESTED_TYPE CGAL_Interval_nt          // For tst4.C

#define CGAL_IA_NO_WARNINGS
#define CGAL_NO_PRECONDITIONS
#include <CGAL/Timer.h>
#include <CGAL/Interval_arithmetic.h>

typedef TESTED_TYPE IA;

// Some simple operators benchmarks.

void bench()
{
  const int loops = 10000000;
  CGAL_Timer t;
  double dt;
  const double dd = 1.0;
  const IA a(0.12);
  // const IA b(2.1);
  const IA b(IA(21)/10);
  IA c(1);

// egcs-1.1 + -O + not advanced n'affiche pas pareil....
// De même que le snapshot du 19 octobre...  inquiétant.
// Quand j'aurai du temps: faire un bug report.
   c = a + b;
   cout << a << b << endl;
   cout << c <<  endl;

  cout << loops << " loops.\n";

#define BENCH_MACRO(op) { \
  dt = t.time(); t.start(); c = 1; \
  for (int i=0; i<loops; i++) { c = a op b; } \
  t.stop(); \
  cout << c << "\t" #op " " << t.time()-dt << endl; \
}

#define BENCH_MACRO_eq(op1,op2) { \
  dt = t.time(); t.start(); c = 1; \
  for (int i=0; i<loops; i++) { c op1 b; c op2 b; } \
  t.stop(); \
  cout << c << "\t" #op1 " " #op2 " " << t.time()-dt << endl; \
}

  BENCH_MACRO (+);
  BENCH_MACRO (*);
  BENCH_MACRO (/);

  BENCH_MACRO_eq (+=, +=);
  BENCH_MACRO_eq (+=, -=);
  BENCH_MACRO_eq (*=, /=);

  dt = t.time(); t.start();
  for (int i=0; i<loops; i++) { c = sqrt(b); }
  t.stop();
  cout << c << "\tsqrt\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (int i=0; i<loops; i++) { c = b; }
  t.stop();
  cout << c << "\t=\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (int i=0; i<loops; i++) { c = c * dd; }
  t.stop();
  cout << c << "\tia*dbl\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (int i=0; i<loops; i++) { c = dd * c; }
  t.stop();
  cout << c << "\tdbl*ia\t" << t.time()-dt << endl;
}


int main()
{
#ifdef ADVANCED
  CGAL_FPU_set_rounding_to_infinity();
  cout << "Benching the class CGAL_Interval_nt_advanced.\n";
#else
  cout << "Benching the class CGAL_Interval_nt.\n";
#endif

  cout.precision(20);
  bench();

#ifdef ADVANCED
  CGAL_FPU_set_rounding_to_nearest();
#endif

  return 0;
}
