// Generic test file for the IA package. include/tst_generic.C
//  $Revision$
//  $Date$
// Written by Sylvain Pion, 1997-1999.

// This file is included from tst[12].C, that do just a #define:
// #define TESTED_TYPE CGAL::Interval_nt_advanced // For tst1.C
// #define TESTED_TYPE CGAL::Interval_nt          // For tst2.C

#define CGAL_IA_NO_EXCEPTION
#define CGAL_IA_NO_WARNINGS
#include <CGAL/Interval_arithmetic.h>

// #define DEBUG(a) a;
#define DEBUG(a)

typedef TESTED_TYPE IA;

// This test program computes the coordinates of a sequence of points drawing
// a spiral.  It tests, using Interval Arithmetic, whether we fall back on an
// axis.  With double precision, the first possible solution is 396.

int spiral_test()
{
  int i=0;
  IA x_i (1), y_i (0), x_ip1, y_ip1, length;

  // try {
  while (++i < 500)
  {
    x_ip1 = x_i - y_i/sqrt((IA)i);
    y_ip1 = y_i + x_i/sqrt((IA)i);
    x_i = x_ip1;
    y_i = y_ip1;
    length = x_i*x_i + y_i*y_i;
    DEBUG(cout<<i<<": (" << x_i << " , " << y_i << ") : " << length << "\n";)
    // if ((x_i == 0) || (y_i == 0))
    if ( x_i.overlap(0) || y_i.overlap(0) )
      break;
  };
  // }
  // catch (IA::unsafe_comparison) { }

  return (i == 396);
}

// Here we iteratively compute sqrt(interval), where interval is [0.5;1.5]
// at the beginning.  It must converge to the fixed point [1-2^-52 ; 1+2^-52].
// NB: Note that 1-2^-52 is not the closest to 1 (which is 1-2^-53)... Funny.

int square_root_test()
{
  int i=0;
  IA a (0.5, 1.5);

  while (++i < 500)
  {
    IA b = sqrt(a);
    DEBUG( cout << a-1 << endl; )
    if ( b.is_same(a) )
      break;
    a = b;
  };
  DEBUG( cout << i; )
  a -= 1;
  return ( (i==54) &&
           (a.upper_bound() == - a.lower_bound()) &&
           (a.upper_bound() == 1/(1024.0*1024*1024*1024*1024*4)) );
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
  IA a (2), b(2.1);
  IA c (-2,2), d(-2.1,2.1);

  DEBUG( cout << "+infinity = " << HUGE_VAL; )
  DEBUG( cout << "  maxdouble = " << IA::max_double << endl; )
  DEBUG( cout << "CGAL::largest = " << IA::largest << endl; )
  DEBUG( cout << "CGAL::smallest = " << IA::smallest << endl; )
  while (++i < 20)
  {
    a *= a;
    b = b * b;
    c *= c;
    d = d * d;
    DEBUG( cout << a << b << c << d << endl; )
    // DEBUG( cout << a << endl; )
  }

  return a.is_same(IA(IA::max_double, HUGE_VAL)) &&
         b.is_same(IA(IA::max_double, HUGE_VAL)) &&
         c.is_same(IA::largest) &&
         d.is_same(IA::largest);
}


// Here we take a initial interval, and we multiply it by itself...
//     The fixed point must be:
// Same thing for [0.5;0.5]      -> [0;MIN_DOUBLE].
// Same thing for [-0.5;0.5]     -> [-MIN_DOUBLE;MIN_DOUBLE].

int underflow_test()
{
  int i=0;
  IA a (0.5), b(-0.5,0.5);

  DEBUG( cout << IA::min_double << endl;)
  while (++i < 20)
  {
    a *= a;
    b = b * b;
    DEBUG( cout << a << b << endl; )
  }

  return a.is_same(IA(0, IA::min_double)) && b.is_same(IA::smallest);
}


// Here we specifically test the division code.
// We iterate the function f(x)= (1/x + x)/4 + 1/2.

int division_test()
{
  IA a (1), b(0);
  IA c = a/b;
  IA d = IA(-1,1)/-2+1; // aka (0.5,1.5);
  IA e (-d);
  int i=0;

  while (++i < 100)
  {
    b = ((IA)1/d + d)/4 + 0.5;
    a = ((IA)-1/e -e*1)/-4 - 0.5; // make it complicated to test more cases.
    DEBUG( cout << d << e << endl; )
    if ( b.is_same(d) && a.is_same(e) )
      break;
    d = b;
    e = a;
  }
  DEBUG( cout << d << e << i << endl; )
  DEBUG( cout << d-1 << e+1 << endl; )

  return c.is_same(IA::largest) && (i == 54);
}


// Here it's mainly to have a 100% coverage for the test-suite.

int multiplication_test()
{
  const IA a (-2,-1), b (-1,1);
  const IA d (-2,2), e (1,2), f (-2,-1);
  IA c, g, h, i, j;
  c = a * b;
  g = d * e;
  h = d * f;
  i = a * e;
  j = j;

  return 1;
}

// Here we test the specialized functions for IA.
// They are usually templated in CGAL, but I've overriden them.

int utility_test()
{
  bool tmpflag, flag = true;
  const IA a(-1,1), b(-1,0), c(0,0), d(0,1), e(1,2), f(-2,-1), g(1);
  IA h = (IA)1/c;

  tmpflag = (CGAL::sign(c) == CGAL::ZERO) &&
            (CGAL::sign(e) == CGAL::POSITIVE) &&
            (CGAL::sign(f) == CGAL::NEGATIVE) ;
  DEBUG( cout << "Sign test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::abs(a).is_same(IA(0,1)) && CGAL::abs(b).is_same(IA(0,1)) &&
            CGAL::abs(c).is_same(IA(0,0)) && CGAL::abs(d).is_same(IA(0,1)) &&
            CGAL::abs(e).is_same(IA(1,2)) && CGAL::abs(f).is_same(IA(1,2)) &&
            CGAL::abs(g).is_same(g) ;
  DEBUG( cout << "CGAL::abs test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::is_valid(a) && CGAL::is_valid(h);
  DEBUG( cout << "CGAL::is_valid test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::is_finite(a) && !CGAL::is_finite(h);
  DEBUG( cout << "CGAL::is_finite test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::max(a,d).is_same(IA(0,1));
  DEBUG( cout << "CGAL::max test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::min(a,b).is_same(IA(-1,0));
  DEBUG( cout << "CGAL::max test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::sign(f) == CGAL::NEGATIVE;
  DEBUG( cout << "CGAL::sign test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  tmpflag = (CGAL::compare(b,e) == CGAL::SMALLER)
         && (CGAL::compare(g,g) == CGAL::EQUAL);
  DEBUG( cout << "CGAL::compare test :\t" << tmpflag << endl; )
  flag = flag && tmpflag;

  return flag;
}


int main()
{
  // unsigned short rd_mode;
  // GETFPCW(rd_mode);
  // cout << hex << rd_mode << endl;

#ifdef ADVANCED
  CGAL::FPU_set_rounding_to_infinity();
  cout << "Stress-testing the class CGAL::Interval_nt_advanced.\n";
#else
  cout << "Stress-testing the class CGAL::Interval_nt.\n";
#endif

#if 0
// Note: it works with -O3, but not -g... (at least, with latest eg++) !?!?!?
  double a = 2.1; // Check 2.0 too.
  a = a*a; a = a*a; a = a*a; a = a*a; a = a*a; a = a*a; a = a*a; a = a*a;
  double b = -a;
  b = b*a;
  b = b*a;
  b = b*a;
  b = b*a;
  IA c = 2.1;
  c *= c; c *= c; c *= c; c *= c; c *= c; c *= c; c *= c; c *= c;
  c *= c; c *= c; c *= c; c *= c; c *= c; c *= c; c *= c; c *= c;
  IA d = 2.0;
  d += d; d += d; d += d; d += d; d += d; d += d; d += d; d += d;
  d += d; d += d; d += d; d += d; d += d; d += d; d += d; d += d;
  d += d; d += d; d += d; d += d; d += d; d += d; d += d; d += d;
  d += d; d += d; d += d; d += d; d += d; d += d; d += d; d += d;
  d += d; d += d; d += d; d += d; d += d; d += d; d += d; d += d;
  d += d; d += d; d += d; d += d; d += d; d += d; d += d; d += d;
#endif

  bool tmpflag, flag = true;
  cout.precision(20);

#if 0
cout << a << "  " << CGAL::is_finite(a) << CGAL::is_valid(a) << endl;
cout << b << "  " << CGAL::is_finite(b) << CGAL::is_valid(b) << endl;
cout << c << "  " << CGAL::is_finite(c) << CGAL::is_valid(c);
cout << CGAL::is_finite(c.lower_bound()) << CGAL::is_valid(c.lower_bound());
cout << CGAL::is_finite(c.upper_bound()) << CGAL::is_valid(c.upper_bound())<<endl;
cout << d << "  " << CGAL::is_finite(d) << CGAL::is_valid(d);
cout << CGAL::is_finite(d.lower_bound()) << CGAL::is_valid(d.lower_bound());
cout << CGAL::is_finite(d.upper_bound()) << CGAL::is_valid(d.upper_bound())<<endl;
#endif

  cout << "Printing test:" << endl;
  cout << (IA)-.7 << endl << (IA)7/10 << endl << (IA)1/0 << endl;

  cout << "Do square_root_test()   \t";
  tmpflag = square_root_test();
  cout << (int) tmpflag << endl;
  flag = tmpflag && flag;

  // GETFPCW(rd_mode);
  // cout << hex << rd_mode << endl;

  cout << "Do spiral_test()        \t";
  tmpflag = spiral_test();
  cout << (int) tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do overflow_test()      \t";
  tmpflag = overflow_test();
  cout << (int) tmpflag << endl;
#ifndef __i386
  flag = tmpflag && flag;
#endif

  cout << "Do underflow_test()     \t";
  tmpflag = underflow_test();
  cout << (int) tmpflag << endl;
#ifndef __i386
  flag = tmpflag && flag;
#endif

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

  cout << (int) (0.0 < IA(1)) << endl;

#ifdef ADVANCED
  CGAL::FPU_set_rounding_to_nearest();
#endif

  return !flag;
}
