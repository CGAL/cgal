// Generic test file for the IA package. include/tst_generic.C
//  $Revision$
//  $Date$
// Written by Sylvain Pion, 1997/1998.

// This file is included from tst[12].C, that do just a #define:
// #define TESTED_TYPE CGAL_Interval_nt_advanced // For tst1.C
// #define TESTED_TYPE CGAL_Interval_nt          // For tst2.C

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

  while (++i < 500)
  {
    x_ip1 = x_i - y_i/sqrt((IA)i);
    y_ip1 = y_i + x_i/sqrt((IA)i);
    x_i = x_ip1;
    y_i = y_ip1;
    length = x_i*x_i + y_i*y_i;
    DEBUG( cout << i << ": (" << x_i << " , " << y_i << ") : " << length << "\n";)
    if ((x_i == 0) || (y_i == 0))
      break;
  };

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
    if ( (b.lower_bound() == a.lower_bound()) &&
         (b.upper_bound() == a.upper_bound()) )
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

  while (++i < 20)
  {
    a *= a;
    b = b * b;
    c *= c;
    d = d * d;
    DEBUG( cout << a << b << c << d << endl; )
  }

  return ( (a.lower_bound() !=  HUGE_VAL) && (a.upper_bound() == HUGE_VAL) &&
           (b.lower_bound() !=  HUGE_VAL) && (b.upper_bound() == HUGE_VAL) &&
           (c.lower_bound() == -HUGE_VAL) && (c.upper_bound() == HUGE_VAL) &&
           (d.lower_bound() == -HUGE_VAL) && (d.upper_bound() == HUGE_VAL) );
}


// Here we take a initial interval, and we multiply it by itself...
//     The fixed point must be:
// Same thing for [0.5;0.5]      -> [0;MIN_DOUBLE].
// Same thing for [-0.5;0.5]     -> [-MIN_DOUBLE;MIN_DOUBLE].

int underflow_test()
{
  int i=0;
  IA a (0.5), b(-0.5,0.5);

  while (++i < 20)
  {
    a *= a;
    b = b * b;
    DEBUG( cout << a << b << endl; )
  }

  return ( (a.lower_bound() ==  0) && (a.upper_bound() != 0) &&
           (b.lower_bound() !=  0) && (b.upper_bound() != 0) );
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
    if ( (b.lower_bound() == d.lower_bound()) &&
         (b.upper_bound() == d.upper_bound()) &&
         (a.upper_bound() == e.upper_bound()) &&
         (a.upper_bound() == e.upper_bound()) )
      break;
    d = b;
    e = a;
  }
  DEBUG( cout << d << e << i << endl; )
  DEBUG( cout << d-1 << e+1 << endl; )

  return ( (c.lower_bound() == -HUGE_VAL) &&
           (c.upper_bound() ==  HUGE_VAL) && (i == 54) );
}


// Here it's just to have a 100% coverage for the test-suite.

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

  return -1;
}


int main()
{
#ifdef ADVANCED
  CGAL_FPU_set_rounding_to_infinity();
  cout << "Stress-testing the class CGAL_Interval_nt_advanced.\n";
#else
  cout << "Stress-testing the class CGAL_Interval_nt.\n";
#endif

  bool tmpflag, flag = true;
  cout.precision(20);

  cout << "Printing test:" << endl;
  cout << (IA)-.7 << endl << (IA)7/10 << endl << (IA)1/0 << endl;

  cout << "Do square_root_test() ";
  tmpflag = square_root_test();
  cout << tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do spiral_test() ";
  tmpflag = spiral_test();
  cout << tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do overflow_test() ";
  tmpflag = overflow_test();
  cout << tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do underflow_test() ";
  tmpflag = underflow_test();
  cout << tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do division_test() ";
  tmpflag = division_test();
  cout << tmpflag << endl;
  flag = tmpflag && flag;

  cout << "Do multiplication_test() ";
  tmpflag = multiplication_test();
  cout << tmpflag << endl;
  flag = tmpflag && flag;

#ifdef ADVANCED
  CGAL_FPU_set_rounding_to_nearest();
#endif

  return !flag;
}
