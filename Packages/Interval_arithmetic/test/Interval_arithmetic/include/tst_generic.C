// Generic test file for the IA package.
// Sylvain Pion, 1997-2000.

// This file is included from tst[12].C, that do just a #define:
// #define TESTED_TYPE Interval_nt_advanced // For tst1.C
// #define TESTED_TYPE Interval_nt          // For tst2.C

#include <CGAL/Interval_arithmetic.h>

// #define DEBUG(a) a;
#define DEBUG(a)

typedef TESTED_TYPE IA_nt;

// The following wrappers are needed, since we can't call the right min/max
// from outside namespace CGAL directly.
namespace CGAL {
template <class FT>
inline FT my_max (const FT &a, const FT &b)
{
  return max(a,b);
}

template <class FT>
inline FT my_min (const FT &a, const FT &b)
{
  return min(a,b);
}
}                                                                               

double my_sqrt(double d)
{
  return CGAL_BUG_SQRT(d);
}

void empty_handler(const char*, const char*, const char*, int, const char *)
{
  // Do nothing.
}

double test_force_to_double(double d, double e)
{
  return CGAL_IA_FORCE_TO_DOUBLE(d+e);
}

IA_nt test_add(const IA_nt &a, const IA_nt &b)
{
  return a+b;
}

IA_nt test_mul_2_add(const IA_nt &a)
{
  return a+a;
}

IA_nt test_mult_2_hand_optimized (const IA_nt &a)
{
  return IA_nt(-CGAL_IA_FORCE_TO_DOUBLE(2.0*(-a.inf())),
                CGAL_IA_FORCE_TO_DOUBLE(2.0*a.sup()));
}

IA_nt test_mult_2_hand_optimized_bis (const IA_nt &a)
{
  return IA_nt(-CGAL_IA_FORCE_TO_DOUBLE((-2.0)*a.inf()),
                CGAL_IA_FORCE_TO_DOUBLE(2.0*a.sup()));
}

IA_nt test_mult_4_hand_optimized (const IA_nt &a)
{
  return IA_nt(-CGAL_IA_FORCE_TO_DOUBLE(4.0*(-a.inf())),
                CGAL_IA_FORCE_TO_DOUBLE(4.0*a.sup()));
}

IA_nt test_mult_4_hand_optimized_bis (const IA_nt &a)
{
  return IA_nt(-CGAL_IA_FORCE_TO_DOUBLE((-4.0)*a.inf()),
                CGAL_IA_FORCE_TO_DOUBLE(4.0*a.sup()));
}

#if 0
IA_nt test_mult_2a(const IA_nt &a) { return 2.0 * a; }
IA_nt test_mult_a2(const IA_nt &a) { return a * 2.0; }

IA_nt test_mult_4a(const IA_nt &a) { return 4.0 * a; }
IA_nt test_mult_a4(const IA_nt &a) { return a * 4.0; }
#endif // 0

IA_nt test_square(const IA_nt &a) { return CGAL_NTS square(a); }
IA_nt test_sqr(const IA_nt &a) { return a*a; }

// This test program computes the coordinates of a sequence of points drawing
// a spiral.  It tests, using Interval Arithmetic, whether we fall back on an
// axis.  With double precision, the first possible solution is 396.

bool spiral_test()
{
  int i=0;
  IA_nt x_i (1.0), y_i (0.0);

  while (++i < 500)
  {
    IA_nt x_ip1 = x_i - y_i/CGAL::sqrt(IA_nt(i));
    IA_nt y_ip1 = y_i + x_i/CGAL::sqrt(IA_nt(i));
    x_i = x_ip1;
    y_i = y_ip1;
    DEBUG( IA_nt length = CGAL_NTS square(x_i) + CGAL_NTS square(y_i); )
    DEBUG(std::cout<<i<<": (" << x_i << " , " << y_i << ") : " << length << "\n";)
    if ( x_i.do_overlap(0) || y_i.do_overlap(0) )
      break;
  };

  return i == 396;
}

// Here we iteratively compute sqrt(interval), where interval is [0.5;1.5]
// at the beginning.  It must converge to the fixed point [1-2^-52 ; 1+2^-52].
// NB: Note that 1-2^-52 is not the closest to 1 (which is 1-2^-53)... Funny.

bool square_root_test()
{
  int i=0;
  IA_nt a (0.5, 1.5);

  while (++i < 500)
  {
    IA_nt b = CGAL::sqrt(a);
    DEBUG ( std::cout << a-1.0 << std::endl; )
    if ( b.is_same(a) )
      break;
    a = b;
  };
  a -= 1.0;
  DEBUG (
  std::cout << "i          = " << i << std::endl;
  std::cout << "sup = -inf : " << (a.sup() == -a.inf()) << std::endl;
  std::cout << "width ok ? : " << (-a.inf() == 1/(double(1<<30)*(1<<22))) << std::endl;
  ) // DEBUG
  return i==54 && a.sup() == - a.inf() && a.sup() == 1/(double(1<<30)*(1<<22));
}


// Here we take an initial interval, and we multiply it by itself...
//     The fixed point must be:
// Same thing for [2;2]      -> [MAX_DOUBLE;+inf].
// Same thing for [2.1;2.1]  -> [MAX_DOUBLE;+inf].
// Same thing for [-2;2]     -> [-inf;+inf].
// Same thing for [-2.1;2.1] -> [-inf;+inf].

bool overflow_test()
{
  int i;
  IA_nt a (2), b(2.1);
  IA_nt c (-2,2), d(-2.1,2.1);
  IA_nt e (-2,2), f(2), g(-2);

  DEBUG( std::cout << "+infinity = " << HUGE_VAL << std::endl; )
  DEBUG( std::cout << "maxdouble = " << CGAL_IA_MAX_DOUBLE << std::endl; )
  DEBUG( std::cout << "largest   = " << CGAL::Interval_nt_advanced::Largest << std::endl; )
  DEBUG( std::cout << "smallest  = " << CGAL::Interval_nt_advanced::Smallest << std::endl; )
  for (i=0; i<20; i++)
  {
    a *= a;
    b = b * b;
    c *= c;
    d = d * d;
    // DEBUG( std::cout << a << b << c << d << std::endl; )
    DEBUG( std::cout << a << std::endl; )
  }

  for (i=0; i<10000; i++)
  {
    e += e;
    f = f+f;
    g += g;
    DEBUG( std::cout << "f = " << f << std::endl; )
  }

  return a.is_same(IA_nt(CGAL_IA_MAX_DOUBLE, HUGE_VAL)) &&
         b.is_same(IA_nt(CGAL_IA_MAX_DOUBLE, HUGE_VAL)) &&
         c.is_same(IA_nt::Largest) &&
         d.is_same(IA_nt::Largest) &&
	 e.is_same(IA_nt::Largest) &&
	 f.is_same(IA_nt(CGAL_IA_MAX_DOUBLE, HUGE_VAL)) &&
	 g.is_same(-f);
}


// Here we take a initial interval, and we multiply it by itself...
//     The fixed point must be:
// Same thing for [0.5;0.5]      -> [0;MIN_DOUBLE].
// Same thing for [-0.5;0.5]     -> [-MIN_DOUBLE;MIN_DOUBLE].

bool underflow_test()
{
  int i;
  IA_nt a(0.5), b(-0.5,0.5), c(0.5);

  for (i=0; i<20; i++) a *= a;
  for (i=0; i<20; i++) b = b * b;
  for (i=0; i<20; i++) c = CGAL_NTS square(c);

  return a.is_same(IA_nt(0, CGAL_IA_MIN_DOUBLE))
      && b.is_same(IA_nt::Smallest)
      && c.is_same(IA_nt(0, CGAL_IA_MIN_DOUBLE));
}


// Here we specifically test the division code.
// We iterate the function f(x)= (1/x + x)/4 + 1/2.

bool division_test()
{
  IA_nt a (1), b(0);
  IA_nt c = a/b;
  IA_nt d = IA_nt(-1,1)/-2+1; // aka (0.5,1.5);
  IA_nt e (-d);
  int i=0;

  while (++i < 100)
  {
    b = (IA_nt(1)/d + d)/4 + 0.5;
    a = (IA_nt(-1)/e -e*1)/-4 - 0.5; // make it complicated to test more cases.
    DEBUG( std::cout << d << e << std::endl; )
    if ( b.is_same(d) && a.is_same(e) )
      break;
    d = b;
    e = a;
  }
  DEBUG( std::cout << d << e << i << std::endl; )
  DEBUG( std::cout << d-1 << e+1 << std::endl; )

  return c.is_same(IA_nt::Largest) && i == 54;
}


// Here it's mainly to have a 100% coverage for the test-suite.

bool multiplication_test()
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
  IA_nt k = 1;
  i=2;
  h = k+i; h = i+k; h += k; h += i;
  h = k-i; h = i-k; h -= k; h -= i;
  h = k*i; h = i*k; h *= k; h *= i;
  h = k/i; h = i/k; h /= k; h /= i;

  return true;
}

// Here we test the specialized functions for IA.
// They are usually templated in CGAL, but I've overriden them.

bool utility_test()
{
  bool tmpflag, flag = true;
  const IA_nt a(-1,1), b(-1,0), c(0,0), d(0,1), e(1,2), f(-2,-1), g(1);
  IA_nt h = IA_nt(1)/c;

  tmpflag = (CGAL_NTS sign(c) == CGAL::ZERO) &&
            (CGAL_NTS sign(e) == CGAL::POSITIVE) &&
            (CGAL_NTS sign(f) == CGAL::NEGATIVE) ;
  DEBUG( std::cout << "Sign test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL_NTS abs(a).is_same(IA_nt(0,1)) &&
            CGAL_NTS abs(b).is_same(IA_nt(0,1)) &&
            CGAL_NTS abs(c).is_same(IA_nt(0,0)) &&
	    CGAL_NTS abs(d).is_same(IA_nt(0,1)) &&
            CGAL_NTS abs(e).is_same(IA_nt(1,2)) &&
	    CGAL_NTS abs(f).is_same(IA_nt(1,2)) &&
            CGAL_NTS abs(g).is_same(g) ;
  DEBUG( std::cout << "abs test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::is_valid(a) && CGAL::is_valid(h);
  DEBUG( std::cout << "is_valid test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::is_finite(a) && !CGAL::is_finite(h);
  DEBUG( std::cout << "is_finite test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::my_max(a,d).is_same(IA_nt(0,1));
  DEBUG( std::cout << "max test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::my_min(a,b).is_same(IA_nt(-1,0));
  DEBUG( std::cout << "min test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL_NTS sign(f) == CGAL::NEGATIVE;
  DEBUG( std::cout << "sign test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = (CGAL_NTS compare(b,e) == CGAL::SMALLER)
         && (CGAL_NTS compare(g,g) == CGAL::EQUAL);
  DEBUG( std::cout << "compare test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  return flag;
}

// Test the is_valid() function.

double zero = 0.0; // I put it here to avoid compiler warnings.

bool is_valid_test()
{
  CGAL::Failure_behaviour backup = CGAL::set_error_behaviour(CGAL::CONTINUE);
  CGAL::Failure_function prev    = CGAL::set_error_handler(empty_handler);

  bool tmpflag, flag = true;
  const double inf = 1.0/zero;
  const double nan = 0.0 * inf;
  const IA_nt a(nan, nan), b(0,nan), c(nan, 0), d(1,0);
  const IA_nt e(0,1), f(0,0);

  tmpflag = CGAL::is_valid(nan);
  DEBUG( std::cout << std::endl; )
  DEBUG( std::cout << "is_valid( " << nan << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL::is_valid(zero);
  DEBUG( std::cout << "is_valid( " << zero << " ) = " << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::is_valid(inf);
  DEBUG( std::cout << "is_valid( " << inf << " ) = " << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::is_valid(a);
  DEBUG( std::cout << "is_valid( " << a << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL::is_valid(b);
  DEBUG( std::cout << "is_valid( " << b << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL::is_valid(c);
  DEBUG( std::cout << "is_valid( " << c << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL::is_valid(d);
  DEBUG( std::cout << "is_valid( " << d << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL::is_valid(e);
  DEBUG( std::cout << "is_valid( " << e << " ) = " << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::is_valid(f);
  DEBUG( std::cout << "is_valid( " << f << " ) = " << tmpflag << std::endl; )
  flag = flag && tmpflag;

  CGAL::set_error_handler(prev);
  CGAL::set_error_behaviour(backup);
  return flag;
}

// Test the is_finite() function.

bool test_ieee_equality(double d)
{
  return CGAL::is_finite(d);
  // return d == d;
}

bool is_finite_test()
{
  bool tmpflag, flag = true;
  const double inf = 1.0/zero;
  const IA_nt a(inf, inf), b(-inf,inf), c(-inf, 0), d(0,inf);
  const IA_nt e(0,1), f(0,0);

  DEBUG(
  const double nan = inf-inf;
  using CGAL::is_finite;
  std::cout << "Test de is_finite(double)" << std::endl;
  std::cout << "is_finite( " << inf << " ) = " << is_finite(inf) << std::endl;
  std::cout << "is_finite( " << 0.0 << " ) = " << is_finite(0.0) << std::endl;
  std::cout << "is_finite( " << 1.0 << " ) = " << is_finite(1.0) << std::endl;
  std::cout << "is_finite( " << -1.0 << " ) = " << is_finite(-1.0) << std::endl;
  std::cout << "is_finite( " << -inf << " ) = " << is_finite(-inf) << std::endl;
  std::cout << "is_finite( " << nan << " ) = " << is_finite(nan) << std::endl;
  )

  tmpflag = CGAL::is_finite(a);
  DEBUG( std::cout << std::endl; )
  DEBUG( std::cout << "is_finite( " << a << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL::is_finite(b);
  DEBUG( std::cout << "is_finite( " << b << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL::is_finite(c);
  DEBUG( std::cout << "is_finite( " << c << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL::is_finite(d);
  DEBUG( std::cout << "is_finite( " << d << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL::is_finite(e);
  DEBUG( std::cout << "is_finite( " << e << " ) = " << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::is_finite(f);
  DEBUG( std::cout << "is_finite( " << f << " ) = " << tmpflag << std::endl; )
  flag = flag && tmpflag;

  return flag;
}

void print_res (bool res)
{
  std::cout << (res ? "ok" : "ERROR") << std::endl;
}

int main()
{
#ifdef ADVANCED
  CGAL::FPU_CW_t backup = CGAL::FPU_get_and_set_cw(CGAL_FE_UPWARD);
  std::cout << "Stress-testing the class Interval_nt_advanced.\n";
#else
  std::cout << "Stress-testing the class Interval_nt.\n";
#endif

  bool tmpflag, flag = true;
  std::cout.precision(20);

  std::cout << "Printing test:" << std::endl;
  std::cout << (IA_nt)-.7 << std::endl;
  std::cout << (IA_nt)7/10 << std::endl;
  std::cout << (IA_nt)1/0 << std::endl;

#define TEST_MACRO(fn) \
  std::cout << #fn << "\t\t"; \
  tmpflag = fn(); \
  print_res(tmpflag); \
  flag = tmpflag && flag;

  TEST_MACRO(square_root_test);
  TEST_MACRO(spiral_test);
  TEST_MACRO(overflow_test);
  TEST_MACRO(underflow_test);
  TEST_MACRO(division_test);
  TEST_MACRO(multiplication_test);
  TEST_MACRO(utility_test);
  TEST_MACRO(is_valid_test);
  TEST_MACRO(is_finite_test);

  print_res(IA_nt(0) < IA_nt(1));

#ifdef ADVANCED
  CGAL::FPU_set_cw(backup);
#endif

  return !flag;
}
