// Generic bench file for the IA package. include/bench_generic.C
//  $Revision$
//  $Date$
// Written by Sylvain Pion, 1997-1999.

// This file is included from tst[34].C, that do just a #define:
// #define TESTED_TYPE Interval_nt_advanced // For tst3.C
// #define TESTED_TYPE Interval_nt          // For tst4.C

#define CGAL_IA_NO_EXCEPTION
#define CGAL_IA_NO_WARNINGS
// #define CGAL_NO_ASSERTIONS
#include <CGAL/Timer.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/predicates_on_ftC2.h>

using namespace CGAL;

typedef TESTED_TYPE IA_nt;

// Not called, only used to watch at the assembly code produced.

IA_nt add (IA_nt a, IA_nt b)
{
  return a+b;
}

// Some simple operators benchmarks.

void bench()
{
#ifndef LOOPS
  const int loops = 1000;
#else
  const int loops = LOOPS;
#endif

  int i;
  Timer t;
  double dt;
  const double dd = 1.0000001;
  const IA_nt a(0.12);
  // const IA_nt b(2.1);
  const IA_nt b(IA_nt(21)/10);
  IA_nt c(1), d(-5.0/3), e(-6.0/7), f(7.0/9);

   c = a + b;
   cout << a << endl;
   cout << b << endl;
   if (b.inf() == b.sup())
     cout << "error" << endl;
   cout << c << endl;

  cout << loops << " loops.\n";

#define BENCH_MACRO(op) { \
  dt = t.time(); t.start(); c = 1; \
  for (i=0; i<loops; i++) { c = a op b; } \
  t.stop(); \
  cout << c << "\t" #op " " << t.time()-dt << endl; \
}

#define BENCH_MACRO_eq(op1,op2) { \
  dt = t.time(); t.start(); c = 1; \
  for (i=0; i<loops; i++) { c op1 b; c op2 b; } \
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
  for (i=0; i<loops; i++) { c = sqrt(b); }
  t.stop();
  cout << c << "\tsqrt\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (i=0; i<loops; i++) { c = square(b); }
  t.stop();
  cout << c << "\tsquare\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (i=0; i<loops; i++) { c = b; }
  t.stop();
  cout << c << "\t=\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (i=0; i<loops; i++) { c = c * dd; }
  t.stop();
  cout << c << "\tia*dbl\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (i=0; i<loops; i++) { c = dd * c; }
  t.stop();
  cout << c << "\tdbl*ia\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (i=0; i<loops; i++) { c = c + dd; }
  t.stop();
  cout << c << "\tia+dbl\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (i=0; i<loops; i++) { c = dd + c; }
  t.stop();
  cout << c << "\tdbl+ia\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (i=0; i<loops; i++) { c = c - dd; }
  t.stop();
  cout << c << "\tia-dbl\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (i=0; i<loops; i++) { c = dd - c; }
  t.stop();
  cout << c << "\tdbl-ia\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (i=0; i<loops; i++) { c = dd/c; }
  t.stop();
  cout << c << "\td/ia\t" << t.time()-dt << endl;

  dt = t.time(); t.start();
  for (i=0; i<loops; i++) { c = c/dd; }
  t.stop();
  cout << c << "\tia/d\t" << t.time()-dt << endl;

#if 1
  cout << a<<b<<c<<d<<endl;
  Orientation o;
  dt = t.time(); t.start();
  for (i=0; i<loops; i++)
    o = orientationC2(a,b,c,d,e,f);
  t.stop();
  cout << (int)o << "\tori2\t" << t.time()-dt << endl;
#endif
}


int main()
{
#ifdef ADVANCED
  FPU_CW_t backup = FPU_get_cw();
  FPU_set_cw(FPU_cw_up);
  cout << "Benching the class Interval_nt_advanced.\n";
#else
  cout << "Benching the class Interval_nt.\n";
#endif
  double d;
  if ((((int) &d) & 7) != 0)
    cout << "Benchmark might not be meaningful due to bad alignment" << endl;

  cout.precision(20);
  bench();

  IA_nt a=1, b=2;
  (int) sign(a);
  (int) compare(a,b);
#if 0
  // It would be nice if it emitted a warning, because c is not initialized.
  IA_nt c;
  a = c+c;
  cout << c << a << endl;
#endif

#ifdef ADVANCED
  FPU_set_cw(backup);
#endif

  return 0;
}
