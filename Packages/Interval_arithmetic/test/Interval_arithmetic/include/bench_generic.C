// Generic bench file for the IA package.
// Sylvain Pion, 1997-1999.

// This file is included from tst[34].C, that do just a #define:
// #define TESTED_TYPE Interval_nt_advanced // For tst3.C
// #define TESTED_TYPE Interval_nt          // For tst4.C

#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/predicates/kernel_ftC2.h>

#include <cassert>

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
  CGAL::Timer t;
  double dt;
  const double dd = 1.0000001;
  const IA_nt a(0.12);
  // const IA_nt b(2.1);
  const IA_nt b(IA_nt(21)/10);
  IA_nt c(1), d(-5.0/3), e(-6.0/7), f(7.0/9);

   c = a + b;
   std::cout << a << std::endl;
   std::cout << b << std::endl;
   if (b.is_point())
     std::cout << "error due to constant propagation" << std::endl;
   std::cout << c << std::endl;

  std::cout << loops << " loops.\n";

  dt = t.time(); t.start();
  for (i=0; i<loops; i++) { c = b; }
  t.stop();
  std::cout << c << "\t=\t" << t.time()-dt << std::endl;

#define BENCH_MACRO_generic(init, op1, op2) { \
  dt = t.time(); t.start(); init; \
  for (i=0; i<loops; i++) { op1; } \
  t.stop(); \
  std::cout << c << "\t" << op2 << "\t" << t.time()-dt << std::endl; \
  assert( ! c.is_point()); \
}

#define EMPTY do {} while(0)

  BENCH_MACRO_generic(EMPTY, c = a + b, "+");
  BENCH_MACRO_generic(EMPTY, c = a * b, "*");
  BENCH_MACRO_generic(EMPTY, c = a / b, "/");

  BENCH_MACRO_generic(c = 1, c += b; c += b, "+= +=");
  BENCH_MACRO_generic(c = 1, c += b; c -= b, "+= -=");
  BENCH_MACRO_generic(c = 1, c *= b; c /= b, "*= /=");

  BENCH_MACRO_generic(EMPTY, c = sqrt(b), "sqrt");
  BENCH_MACRO_generic(EMPTY, c = square(b), "square");
  BENCH_MACRO_generic(EMPTY, c = c * dd, "ia*d");
  BENCH_MACRO_generic(EMPTY, c = dd * c, "d*ia");
  BENCH_MACRO_generic(EMPTY, c = c + dd, "ia+d");
  BENCH_MACRO_generic(EMPTY, c = dd + c, "d+ia");
  BENCH_MACRO_generic(EMPTY, c = c - dd, "ia-d");
  BENCH_MACRO_generic(EMPTY, c = dd - c, "d-ia");
  BENCH_MACRO_generic(EMPTY, c = dd / c, "d/ia");
  BENCH_MACRO_generic(EMPTY, c = c / dd, "ia/d");

#if 1
  std::cout << a<<b<<c<<d<<std::endl;
  CGAL::Orientation o;
  dt = t.time(); t.start();
  for (i=0; i<loops; i++)
    o = orientationC2(a,b,c,d,e,f);
  t.stop();
  std::cout << (int)o << "\tori2\t" << t.time()-dt << std::endl;
#endif
}


int main()
{
#ifdef ADVANCED
  CGAL::FPU_CW_t backup = CGAL::FPU_get_cw();
  CGAL::FPU_set_cw(CGAL::FPU_cw_up);
  std::cout << "Benching the class Interval_nt_advanced.\n";
#else
  std::cout << "Benching the class Interval_nt.\n";
#endif
  double d;
  if ((((int) &d) & 7) != 0)
    std::cout << "Benchmark might not be meaningful due to bad alignment\n";

  std::cout.precision(20);
  bench();

  IA_nt a=1, b=2;
  (int) sign(a);
  (int) compare(a,b);
#if 0
  // It would be nice if it emitted a warning, because c is not initialized.
  IA_nt c;
  a = c+c;
  std::cout << c << a << std::endl;
#endif

#ifdef ADVANCED
  CGAL::FPU_set_cw(backup);
#endif

  return 0;
}
