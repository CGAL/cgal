// Test file for the Interval_nt<bool> class.
// Sylvain Pion, 1997-2005.

#include <CGAL/config.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/exceptions.h>

// #define DEBUG(a) a;
#define DEBUG(a)


// This test program computes the coordinates of a sequence of points drawing
// a spiral.  It tests, using Interval Arithmetic, whether we fall back on an
// axis.  With double precision, the first possible solution is 396.

template < typename IA_nt >
bool spiral_test()
{
  int i=0;
  IA_nt x_i (1.0), y_i (0.0);

  while (++i < 500)
  {
    IA_nt x_ip1 = x_i - y_i/CGAL_NTS sqrt(IA_nt(i));
    IA_nt y_ip1 = y_i + x_i/CGAL_NTS sqrt(IA_nt(i));
    x_i = x_ip1;
    y_i = y_ip1;
    DEBUG( IA_nt length = CGAL_NTS square(x_i) + CGAL_NTS square(y_i); )
    DEBUG(std::cout<<i<<": (" << x_i << " , " << y_i << ") : " << length << "\n";)
    if ( x_i.do_overlap(0) || y_i.do_overlap(0) )
      break;
  };

#ifdef CGAL_ALWAYS_ROUND_TO_NEAREST
  return i == 365;
#else
  return i == 396;
#endif
}

// Tests for constant propagation through intervals.
// This must not be performed otherwise rounding modes are ignored.
// On the other hand, if we always round to nearest, then constant propagation
// is desirable.
// Note: Non-inlined operators usually stop cprop (*, /, sqrt).
template < typename IA_nt >
bool cprop_test()
{
  // Testing cprop through +.
  IA_nt add = IA_nt(0.00001)+10.1;
  bool good_add = !add.is_point();
  if (!good_add)
    std::cerr << "ERROR : Constant propagation through operator+." <<std::endl;

  // Testing cprop through -.
  IA_nt sub = IA_nt(0.00001)-10.1;
  bool good_sub = !sub.is_point();
  if (!good_sub)
    std::cerr << "ERROR : Constant propagation through operator-." <<std::endl;

  // Testing cprop through *.
  IA_nt mul = IA_nt(0.00001)*10.1;
  bool good_mul = !mul.is_point();
  if (!good_mul)
    std::cerr << "ERROR : Constant propagation through operator*." <<std::endl;

  // Testing cprop through /.
  IA_nt div = IA_nt(0.00001)/10.1;
  bool good_div = !div.is_point();
  if (!good_div)
    std::cerr << "ERROR : Constant propagation through operator/." <<std::endl;

  // Testing cprop through sqrt.
  IA_nt sqrt2 = CGAL_NTS sqrt(IA_nt(2));
  bool good_sqrt = !sqrt2.is_point();
  if (!good_sqrt)
    std::cerr << "ERROR : Constant propagation through sqrt()." <<std::endl;

  return good_add && good_sub && good_mul && good_div && good_sqrt;
}

// Here we iteratively compute sqrt(interval), where interval is [0.5;1.5]
// at the beginning.  It must converge to the fixed point [1-2^-52 ; 1+2^-52].
// NB: Note that 1-2^-52 is not the closest to 1 (which is 1-2^-53)... Funny.

template < typename IA_nt >
bool square_root_test()
{
  int i=0;
  IA_nt a (0.5, 1.5);

  while (++i < 500)
  {
    IA_nt b = CGAL_NTS sqrt(a);
    DEBUG ( std::cout << a-1.0 << std::endl; )
    if ( b.is_same(a) )
      break;
    a = b;
  };
  a -= 1.0;
#ifdef CGAL_ALWAYS_ROUND_TO_NEAREST
  DEBUG (
  std::cout << "i          = " << i << std::endl;
  std::cout << "sup        : " << a.sup() << std::endl;
  std::cout << "inf        : " << a.inf() << std::endl;
  ) // DEBUG
  if (i != 54) {
    return false;
  }
  // When we round to nearest it doesn't quite converge.
  if (a.sup() > 3/(double(1<<30)*(1<<22))) {
    return false;
  }
  if (-3/(double(1<<30)*(1<<22)) > a.inf()) {
    return false;
  }
  return true;
#else
  DEBUG (
  std::cout << "i          = " << i << std::endl;
  std::cout << "sup = -inf : " << (a.sup() == -a.inf()) << std::endl;
  std::cout << "width ok ? : " << (-a.inf() == 1/(double(1<<30)*(1<<22))) << std::endl;
  ) // DEBUG
  return i==54 && a.sup() == - a.inf() && a.sup() == 1/(double(1<<30)*(1<<22));
#endif
}


// Here we take an initial interval, and we multiply it by itself...
//     The fixed point must be:
// Same thing for [2;2]      -> [MAX_DOUBLE;+inf].
// Same thing for [2.1;2.1]  -> [MAX_DOUBLE;+inf].
// Same thing for [-2;2]     -> [-inf;+inf].
// Same thing for [-2.1;2.1] -> [-inf;+inf].

template < typename IA_nt >
bool overflow_test()
{
  int i;
  IA_nt a (2), b(2.1);
  IA_nt c (-2,2), d(-2.1,2.1);
  IA_nt e (-2,2), f(2), g(-2);

  DEBUG( std::cout << "+infinity = " << std::numeric_limits<double>::infinity() << std::endl; )
  DEBUG( std::cout << "maxdouble = " << CGAL_IA_MAX_DOUBLE << std::endl; )
  DEBUG( std::cout << "largest   = " << CGAL::Interval_nt_advanced::largest() << std::endl; )
  DEBUG( std::cout << "smallest  = " << CGAL::Interval_nt_advanced::smallest() << std::endl; )

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

  return a.is_same(IA_nt(CGAL_IA_MAX_DOUBLE, std::numeric_limits<double>::infinity())) &&
         b.is_same(IA_nt(CGAL_IA_MAX_DOUBLE, std::numeric_limits<double>::infinity())) &&
         c.is_same(IA_nt::largest()) &&
         d.is_same(IA_nt::largest()) &&
         e.is_same(IA_nt::largest()) &&
         f.is_same(IA_nt(CGAL_IA_MAX_DOUBLE, std::numeric_limits<double>::infinity())) &&
         g.is_same(-f);
}


// Here we take a initial interval, and we multiply it by itself...
//     The fixed point must be:
// Same thing for [0.5;0.5]      -> [0;MIN_DOUBLE].
// Same thing for [-0.5;0.5]     -> [-MIN_DOUBLE;MIN_DOUBLE].

template < typename IA_nt >
bool underflow_test()
{
  int i;
  IA_nt a(0.5), b(-0.5,0.5), c(0.5);

  for (i=0; i<20; i++) a *= a;
  for (i=0; i<20; i++) b = b * b;
  for (i=0; i<20; i++) c = CGAL_NTS square(c);

#ifdef CGAL_ALWAYS_ROUND_TO_NEAREST
  return a.is_same(IA_nt(-CGAL_IA_MIN_DOUBLE, CGAL_IA_MIN_DOUBLE))
      && b.is_same(IA_nt::smallest())
      && c.is_same(IA_nt(0, CGAL_IA_MIN_DOUBLE));
#else
  return a.is_same(IA_nt(0, CGAL_IA_MIN_DOUBLE))
      && b.is_same(IA_nt::smallest())
      && c.is_same(IA_nt(0, CGAL_IA_MIN_DOUBLE));
#endif
}


// Here we specifically test the division code.
// We iterate the function f(x)= (1/x + x)/4 + 1/2.

template < typename IA_nt >
bool division_test()
{
  IA_nt a (1), b(0);
  IA_nt c = a/b;
  IA_nt d = IA_nt(-1,1)/-2+1L; // aka (0.5,1.5);
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

  return c.is_same(IA_nt::largest()) && i == 54;
}


// Here it's mainly to have a 100% coverage for the test-suite.

template < typename IA_nt >
bool multiplication_test()
{
  const IA_nt a (-2,-1), b (-1,1);
  const IA_nt d (-2,2), e (1,2), f (-2,-1);
  IA_nt c, g, h, i, j;
  c = a * b;
  g = d * e;
  h = d * f;
  i = a * e;
  j = (IA_nt&)j;

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
// They are usually templated in CGAL, but I've overridden them.

template < typename IA_nt >
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

  tmpflag = CGAL_NTS is_finite(a) && !CGAL_NTS is_finite(h);
  DEBUG( std::cout << "is_finite test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = (CGAL::max)(a,d).is_same(IA_nt(0,1));
  DEBUG( std::cout << "max test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = (CGAL::min)(a,b).is_same(IA_nt(-1,0));
  DEBUG( std::cout << "min test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL_NTS sign(f) == CGAL::NEGATIVE;
  DEBUG( std::cout << "sign test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = ((CGAL_NTS compare(b,e)) == CGAL::SMALLER)
         && ((CGAL_NTS compare(g,g)) == CGAL::EQUAL);
  DEBUG( std::cout << "compare test :\t" << tmpflag << std::endl; )
  flag = flag && tmpflag;

  return flag;
}

// Test the is_valid() function.

double zero = 0.0; // I put it here to avoid compiler warnings.

template < typename IA_nt >
bool is_valid_test()
{
  bool tmpflag, flag = true;
  const double inf = 1.0/zero;
  const double nan = zero * inf;
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

  tmpflag = CGAL::is_valid(e);
  DEBUG( std::cout << "is_valid( " << e << " ) = " << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL::is_valid(f);
  DEBUG( std::cout << "is_valid( " << f << " ) = " << tmpflag << std::endl; )
  flag = flag && tmpflag;


  // Now do some tests which can trigger assertion failures, hence the try/catch.
  std::cout << "Testing some assertion protected codes (possible error messages)" << std::endl;
  try {
    const IA_nt a(nan, nan), b(0,nan), c(nan, 0), d(1,0);

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
  }
  catch (CGAL::Assertion_exception&)
  {}

  return flag;
}

// Test the is_finite() function.

template < typename IA_nt >
bool is_finite_test()
{
  bool tmpflag, flag = true;
  const double inf = 1.0/zero;
  const IA_nt a(inf, inf), b(-inf,inf), c(-inf, 0), d(0,inf);
  const IA_nt e(0,1), f(0,0);

  DEBUG(
  const double nan = inf-inf;
  using CGAL_NTS is_finite;
  std::cout << "Test de is_finite(double)" << std::endl;
  std::cout << "is_finite( " << inf << " ) = " << is_finite(inf) << std::endl;
  std::cout << "is_finite( " << 0.0 << " ) = " << is_finite(0.0) << std::endl;
  std::cout << "is_finite( " << 1.0 << " ) = " << is_finite(1.0) << std::endl;
  std::cout << "is_finite( " << -1.0 << " ) = " << is_finite(-1.0) << std::endl;
  std::cout << "is_finite( " << -inf << " ) = " << is_finite(-inf) << std::endl;
  std::cout << "is_finite( " << nan << " ) = " << is_finite(nan) << std::endl;
  )

  tmpflag = CGAL_NTS is_finite(a);
  DEBUG( std::cout << std::endl; )
  DEBUG( std::cout << "is_finite( " << a << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL_NTS is_finite(b);
  DEBUG( std::cout << "is_finite( " << b << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL_NTS is_finite(c);
  DEBUG( std::cout << "is_finite( " << c << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL_NTS is_finite(d);
  DEBUG( std::cout << "is_finite( " << d << " ) = " << tmpflag << std::endl; )
  flag = flag && !tmpflag;

  tmpflag = CGAL_NTS is_finite(e);
  DEBUG( std::cout << "is_finite( " << e << " ) = " << tmpflag << std::endl; )
  flag = flag && tmpflag;

  tmpflag = CGAL_NTS is_finite(f);
  DEBUG( std::cout << "is_finite( " << f << " ) = " << tmpflag << std::endl; )
  flag = flag && tmpflag;

  return flag;
}

void print_res (bool res)
{
  std::cout << (res ? "ok" : "ERROR") << std::endl;
}

template < typename IA_nt >
bool test ()
{
  typename IA_nt::Protector p;

  bool tmpflag, flag = true;
  std::cout.precision(20);

  std::cout << "Printing test:" << std::endl;
  std::cout << (IA_nt)-.7 << std::endl;
  std::cout << (IA_nt)7/10 << std::endl;
  std::cout << (IA_nt)1/0 << std::endl;

#define TEST_MACRO(fn) \
  std::cout << #fn << "\t\t"; \
  tmpflag = fn<IA_nt>(); \
  print_res(tmpflag); \
  flag = tmpflag && flag;

  TEST_MACRO(cprop_test);
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

  return flag;
}

int main()
{
#ifdef CGAL_ALWAYS_ROUND_TO_NEAREST
  std::cout << "Stress-testing the class Interval_nt<> always rounding to nearest.\n";
  bool ok = test<CGAL::Interval_nt<> >();
  std::cout << "\nStress-testing the class Interval_nt_advanced always rounding to nearest.\n";
  ok &= test<CGAL::Interval_nt_advanced>();
#else
  std::cout << "Stress-testing the class Interval_nt<>.\n";
  bool ok = test<CGAL::Interval_nt<> >();
  std::cout << "\nStress-testing the class Interval_nt_advanced.\n";
  ok &= test<CGAL::Interval_nt_advanced>();
#endif

  return !ok;
}
