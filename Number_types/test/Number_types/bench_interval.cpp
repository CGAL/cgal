// Generic bench file for the IA package.
// Sylvain Pion, 1997-2005.

#include <CGAL/Timer.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/predicates/kernel_ftC2.h>

#include <CGAL/disable_warnings.h>

#include <cassert>

#ifndef LOOPS
#  define LOOPS 1000
#endif

const int loops = LOOPS;

// Some simple operators benchmarks.

template < typename IA_nt >
void bench()
{
  int i;
  CGAL::Timer t;
  double dt;
  const double dd = 1.0000001;
  const IA_nt a(0.12);
  // const IA_nt b(2.1);
  IA_nt b(IA_nt(21.0)/10.0);
  IA_nt c(1);

  c = a + b;
  std::cout << a << std::endl << b << std::endl;
  if (b.is_point())
  {
    std::cout << "error due to constant propagation" << std::endl;
    assert(false);
  }
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

  BENCH_MACRO_generic(EMPTY, c = CGAL_NTS sqrt(b), "sqrt");
  BENCH_MACRO_generic(EMPTY, c = CGAL_NTS square(b), "square");
  BENCH_MACRO_generic(EMPTY, c = c * IA_nt(dd), "ia*d");
  BENCH_MACRO_generic(EMPTY, c = IA_nt(dd) * c, "d*ia");
  BENCH_MACRO_generic(EMPTY, c = c + IA_nt(dd), "ia+d");
  BENCH_MACRO_generic(EMPTY, c = IA_nt(dd) + c, "d+ia");
  BENCH_MACRO_generic(EMPTY, c = c - IA_nt(dd), "ia-d");
  BENCH_MACRO_generic(EMPTY, c = IA_nt(dd) - c, "d-ia");
  BENCH_MACRO_generic(EMPTY, c = IA_nt(dd) / c, "d/ia");
  BENCH_MACRO_generic(EMPTY, c = c / IA_nt(dd), "ia/d");
}

// OrientationC2() benchmark.

template < typename IA_nt >
void bench_orientation()
{
  IA_nt a(0.12);
  IA_nt b(IA_nt(21.0)/10.0);
  IA_nt c(1), d(-5.0/3), e(-6.0/7), f(7.0/9);

  int i;
  CGAL::Timer t;
  double dt;
  std::cout << a << std::endl << b << std::endl;
  std::cout << c << std::endl << d << std::endl;
  CGAL::Orientation o;
  dt = t.time(); t.start();
  for (i=0; i<loops; i++)
    o = CGAL::orientationC2(a,b,c,d,e,f);
  t.stop();
  std::cout << (int)o << "\tori2\t" << t.time()-dt << std::endl;
}


template < typename IA_nt >
void bench_main()
{
  typename IA_nt::Protector p;

  double d;
  if ((((std::size_t) &d) & 7) != 0)
    std::cout << "Benchmark might not be meaningful due to bad alignment\n";

  std::cout.precision(20);
  bench<IA_nt>();
  bench_orientation<IA_nt>();

  IA_nt a=1, b=2;
  (void) CGAL_NTS sign(a);
  (void) CGAL_NTS compare(a,b);
}

int main()
{
  std::cout << "Benching the class Interval_nt.\n";
  bench_main<CGAL::Interval_nt<> >();
  std::cout << "\nBenching the class Interval_nt_advanced.\n";
  bench_main<CGAL::Interval_nt_advanced>();
  return 0;
}
