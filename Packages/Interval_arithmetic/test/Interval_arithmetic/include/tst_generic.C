// Generic test file for the IA package. include/tst_generic.C
//  $Revision$
//  $Date$
// Written by Sylvain Pion, 1997-1999.

// This file is included from tst[12].C, that do just a #define:
// #define TESTED_TYPE Interval_nt_advanced // For tst1.C
// #define TESTED_TYPE Interval_nt          // For tst2.C

#define CGAL_IA_NO_EXCEPTION
#define CGAL_IA_NO_WARNINGS
#include <CGAL/Interval_arithmetic.h>

using namespace CGAL;

// #define DEBUG(a) a;
#define DEBUG(a)

typedef TESTED_TYPE IA_nt;

// This test program computes the coordinates of a sequence of points drawing
// a spiral.  It tests, using Interval Arithmetic, whether we fall back on an
// axis.  With double precision, the first possible solution is 396.

int spiral_test()
{
  int i=0;
  IA_nt x_i (1), y_i (0), x_ip1, y_ip1, length;

  // try {
  while (++i < 500)
  {
    x_ip1 = x_i - y_i/sqrt((IA_nt)i);
    y_ip1 = y_i + x_i/sqrt((IA_nt)i);
    x_i = x_ip1;
    y_i = y_ip1;
    length = x_i*x_i + y_i*y_i;
    DEBUG(cout<<i<<": (" << x_i << " , " << y_i << ") : " << length << "\n";)
    // if ((x_i == 0) || (y_i == 0))
    if ( x_i.overlap(0) || y_i.overlap(0) )
      break;
  };
  // }
  // catch (IA_nt::unsafe_comparison) { }

  return (i == 396);
}

// Here we iteratively compute sqrt(interval), where interval is [0.5;1.5]
// at the beginning.  It must converge to the fixed point [1-2^-52 ; 1+2^-52].
// NB: Note that 1-2^-52 is not the closest to 1 (which is 1-2^-53)... Funny.

int square_root_test()
{
  int i=0;
  IA_nt a (0.5, 1.5);

  while (++i < 500)
  {
    IA_nt b = sqrt(a);
    DEBUG( cout << a-1 << endl; )
    if ( b.is_same(a) )
      break;
    a = b;
  };
  DEBUG( cout << i; )
  a -= 1;
  return i==54
      && a.sup() == - a.inf()
      && a.sup() == 1/(1024.0*1024*1024*1024*1024*4);
}


// Here we take an initial interval, and we multiply it by itself...
//     The fixed point must be:
// Same thing for [2;2]      -> [MAX_DOUBLE;+inf].
// Same thing for [2.1;2.1]  -> [MAX_DOUBLE;+inf].
// Same thing for [-2;2]     -> [-inf;+inf].
// Same thing for [-2.1;2.1] -> [-inf;+inf].

int overflow_test()
{
  int i=0;
  IA_nt a (2), b(2.1);
  IA_nt c (-2,2), d(-2.1,2.1);
  IA_nt e (-2,2), f(2), g(-2);

  DEBUG( cout << "+infinity = " << HUGE_VAL; )
  DEBUG( cout << "  maxdouble = " << CGAL_IA_MAX_DOUBLE << endl; )
  DEBUG( cout << "largest = " << CGAL_IA_LARGEST << endl; )
  DEBUG( cout << "smallest = " << CGAL_IA_SMALLEST << endl; )
  for (i=0; i<20; i++)
  {
    a *= a;
    b = b * b;
    c *= c;
    d = d * d;
    DEBUG( cout << a << b << c << d << endl; )
    // DEBUG( cout << a << endl; )
  }

  for (i=0; i<10000; i++)
  {
    e += e;
    f = f+f;
    g += g;
  }

  return a.is_same(IA_nt(CGAL_IA_MAX_DOUBLE, HUGE_VAL)) &&
         b.is_same(IA_nt(CGAL_IA_MAX_DOUBLE, HUGE_VAL)) &&
         c.is_same(CGAL_IA_LARGEST) &&
         d.is_same(CGAL_IA_LARGEST) &&
	 e.is_same(CGAL_IA_LARGEST) &&
	 f.is_same(IA_nt(CGAL_IA_MAX_DOUBLE, HUGE_VAL)) &&
	 g.is_same(-f);
}


// Here we take a initial interval, and we multiply it by itself...
//     The fixed point must be:
// Same thing for [0.5;0.5]      -> [0;MIN_DOUBLE].
// Same thing for [-0.5;0.5]     -> [-MIN_DOUBLE;MIN_DOUBLE].

int underflow_test()
{
  int i;
  IA_nt a(0.5), b(-0.5,0.5), c(0.5);

  for (i=0; i<20; i++) a *= a;
  for (i=0; i<20; i++) b = b * b;
  for (i=0; i<20; i++) c = square(c);

  return a.is_same(IA_nt(0, CGAL_IA_MIN_DOUBLE))
      && b.is_same(CGAL_IA_SMALLEST)
      && c.is_same(IA_nt(0, CGAL_IA_MIN_DOUBLE));
}


// Here we specifically test the division code.
// We iterate the function f(x)= (1/x + x)/4 + 1/2.

int division_test()
{
  IA_nt a (1), b(0);
  IA_nt c = a/b;
  IA_nt d = IA_nt(-1,1)/-2+1; // aka (0.5,1.5);
  IA_nt e (-d);
  int i=0;

  while (++i < 100)
  {
    b = (1/d + d)/4 + 0.5;
    a = (-1/e -e*1)/-4 - 0.5; // make it complicated to test more cases.
    DEBUG( cout << d << e << endl; )
    if ( b.is_same(d) && a.is_same(e) )
      break;
    d = b;
    e = a;
  }
  DEBUG( cout << d << e << i << endl; )
  DEBUG( cout << d-1 << e+1 << endl; )

  return c.is_same(CGAL_IA_LARGEST) && (i == 54);
}


// Here it's mainly to have a 100% coverage for the test-suite.

int multiplication_test()
{
  const IA_nt a (-2,-1), b (-1,1);
  const IA_nt d (-2,2), e (1,2), f (-2,-1);
  IA_nt c, g, h, i, j;
  c = a * b;
  g = d * e;
  h = d * f;
  i = a * e;
  j = j;

  // When CGAL_IA_DEBUG is defined, it'll test the current rounding mode for
  // these operations.
  double k = 1;
  i=2;
  h = k+i; h = i+k; h += k; h += i;
  h = k-i; h = i-k; h -= k; h -= i;
  h = k*i; h = i*k; h *= k; h *= i;
  h = k/i; h = i/k; h /= k; h /= i;

  return 1;
}

// Here we test the specialized functions for IA.
// They are usually templated in CGAL, but I've overriden them.

int utility_test()
{
  bool tmpflag, flag = true;
  const IA_nt a(-1,1), b(-1,0), c(0,0), d(0,1), e(1,2), f(-2,-1), g(1);
  IA_nt h = 1/c;

  tmpflag = (sign(c) == ZERO) &&
            (sign(e) == POSITIVE) &&
            (sign(f) == NEGATIVE) ;
  DEBUG( cout << "Sign test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = abs(a).is_same(IA_nt(0,1)) && abs(b).is_same(IA_nt(0,1)) &&
            abs(c).is_same(IA_nt(0,0)) && abs(d).is_same(IA_nt(0,1)) &&
            abs(e).is_same(IA_nt(1,2)) && abs(f).is_same(IA_nt(1,2)) &&
            abs(g).is_same(g) ;
  DEBUG( cout << "abs test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = is_valid(a) && is_valid(h);
  DEBUG( cout << "is_valid test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = is_finite(a) && !is_finite(h);
  DEBUG( cout << "is_finite test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = max(a,d).is_same(IA_nt(0,1));
  DEBUG( cout << "max test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = min(a,b).is_same(IA_nt(-1,0));
  DEBUG( cout << "max test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = sign(f) == NEGATIVE;
  DEBUG( cout << "sign test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = (compare(b,e) == SMALLER)
         && (compare(g,g) == EQUAL);
  DEBUG( cout << "compare test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  return flag;
}


int main()
{
#ifdef ADVANCED
  FPU_CW_t backup = FPU_get_cw();
  FPU_set_cw(FPU_cw_up);
  cout << "Stress-testing the class Interval_nt_advanced.\n";
#else
  cout << "Stress-testing the class Interval_nt.\n";
#endif

  bool tmpflag, flag = true;
  cout.precision(20);

  cout << "Printing test:" << endl;
  cout << (IA_nt)-.7 << endl << (IA_nt)7/10 << endl << (IA_nt)1/0 << endl;

  cout << "Do square_root_test()   \t";
  tmpflag = square_root_test();
  cout << (int) tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do spiral_test()        \t";
  tmpflag = spiral_test();
  cout << (int) tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do overflow_test()      \t";
  tmpflag = overflow_test();
  cout << (int) tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do underflow_test()     \t";
  tmpflag = underflow_test();
  cout << (int) tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do division_test()      \t";
  tmpflag = division_test();
  cout << (int) tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do multiplication_test()\t";
  tmpflag = multiplication_test();
  cout << (int) tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do utility_test()       \t";
  tmpflag = utility_test();
  cout << (int) tmpflag << endl;
  flag = tmpflag && flag;

  cout << (int) (0.0 < IA_nt(1)) << endl;

#ifdef ADVANCED
  FPU_set_cw(backup);
#endif

  return !flag;
}
