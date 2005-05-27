// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Athanasios Kakargias <grad0460@di.uoa.gr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/Root_of/Root_of_2.h

#ifndef CGAL_ROOT_OF_ROOT_OF_2_H
#define CGAL_ROOT_OF_ROOT_OF_2_H

#include <iostream>
#include <cassert>
#include <CGAL/enum.h>
#include <CGAL/NT_extensions_Root_of/Root_of_utils.h>
//#include <CGAL/Root_of/Root_of_1.h>
#include <CGAL/Algebraic_kernel/internal_functions_comparison_root_of_2.h>
#include <CGAL/Interval_arithmetic.h>

#include <cassert>
#include <CGAL/Quotient.h>

namespace CGAL {

// Number type representing a real root of a polynomial
// of degree 1 or 2 over RT.
//
// It supports :
// - constructor from degree 2 polynomial coefficients and a boolean
// - constructor from degree 1 polynomial coefficients
// - constructor from RT
// - unary operator-()
// - additions, subtractions, multiplications with an RT.
// - additions, subtractions, multiplications with an RootOf_1.
// - square()
// - <, >, <=, >=, ==, != (symetric, mixed with RT, mixed with RootOf_1, mixed with int)
// - compare()            (symetric, mixed with RT, mixed with RootOf_1, mixed with int)
// - sign()
// - to_double()
// - to_interval()
// - is_valid()
// - operator[] to access the coefficients  (leading coeff is always positive)
// - .conjuguate()
// - .discriminant()
// - .eval_at()
// - .sign_at()
// - .degree()
// - .is_valid()
// - io_tag()
// - operator<<()
// - print() ("pretty" printing)
// - make_root_of_2()
//
// TODO :
// - Find a number type concept which is independent on namespace CGAL.
//   We already partially do this by only requiring the return type of
//   sign() and compare() to be convertible and comparable to int
//   (and not force it to be some particular enum).
// - replace the duplicated operators by boost's operators.hpp
//   (we should really start planning Boost's integration with CGAL)
// - add inverse(), and division of/by an RT.
// - add subtraction/addition with a degree 2 Root_of of the same field ?
// - add sqrt() (when it's degree 1), or a make_sqrt<RT>(const RT &r) ?
// - add +, -, *, / (when it is possible) ?
// - add constructor from CGAL::Polynomial ?
//   There should be a proper separate class Polynomial.
// - in compare_roots, we evaluate the polynomial at some FT, or at some
//   root of degree 1 polynomials.  It would be nice to have a separate
//   polynomial class which performed this task (and others)...
// - overloaded versions of make_root_of_2<>() for Lazy_exact_nt<> and others.


template < typename RT_ >
class Root_of_2 {
  RT_  C0,C1,C2; // Coefficients (see below)
  char _smaller; // Is it the smaller of the two roots (for degree 2) ?
                 // (we use a char because it's smaller than a bool)

  // the value is the root of P(X) = C2.X^2 + C1.X + C0,
  // and C2 > 0.
  // _smaller indicates if it's the smaller of the 2 roots.

public:

  typedef RT_    RT;

  Root_of_2()
    : C0(0), C1(0), C2(1)
  {
    assert(is_valid());
  }

  Root_of_2(const RT& c0)
    : C0(- square(c0)), C1(0), C2(1), _smaller(c0<0)
  {
    assert(is_valid());
  }

  Root_of_2(const RT& c1, const RT& c0)
    : C0( - square(c0)), C1(0), C2(square(c1)),
      _smaller( sign(c1) == sign(c0) )
  {
    assert(is_valid());
  }

  Root_of_2(const RT& a, const RT& b, const RT& c, const bool s)
    : C0(c), C1(b), C2(a), _smaller(s)
  {
    assert(a != 0);
    if (a < 0) {
      C0 = -c;
      C1 = -b;
      C2 = -a;
    }
    assert(is_valid());
  }

  Root_of_2 operator-() const
  {
    return Root_of_2 (C2, -C1, C0, !_smaller);
  }

  bool is_valid() const
  {
    return (C2 > 0) && ! is_negative(discriminant());
  }

  // The following functions deal with the internal polynomial.
  // Probably they should move to a separate polynomial class.

  bool is_smaller() const
  {
    return _smaller;
  }

  const RT & operator[](int i) const
  {
    assert(i<3);
    return (&C0)[i];
  }

  RT discriminant() const
  {
    return square(C1) - 4*C0*C2;
  }

  template < typename T > 
  T eval_at(const T& x) const 
  { 
    return C0 + x * (C1 + x * C2); 
  } 

  template < typename T > 
  Sign sign_at(const T &x) const
  {
    // Maybe there is slightly more efficient.
    return sign(eval_at(x));
  }

  Root_of_2 conjugate() const
  {
    return Root_of_2(C2, C1, C0, !smaller);
  }
}; // Root_of_2


    template < typename RT >  class Root_of_2;
    template < typename RT >  class Root_of_3;
    template < typename RT >  class Root_of_4;

    // Template default version generating a Root_of_2<>.
    template < typename RT >
    inline
    Root_of_2<RT>
    make_root_of_2(const RT &a, const RT &b, const RT &c, bool smaller)
    {
      assert( a != 0 );
      return Root_of_2<RT>(a, b, c, smaller);
    }

    //Default Traits class for RT types
    template < typename RT >
    struct Root_of_traits
    {
      typedef Quotient< RT >   RootOf_1;
      typedef Root_of_2< RT >  RootOf_2;
      typedef Root_of_3< RT >  RootOf_3;
      typedef Root_of_4< RT >  RootOf_4;
    };

  namespace CGALi {

    // This version is internal and can be re-used for
    // number types which also support division and sqrt().
    template < typename NT >
    NT
    make_root_of_2_sqrt(const NT &a, const NT &b, const NT &c, bool smaller)
    {
      assert( a != 0 );
      NT discriminant = CGAL::square(b) - a*c*4;
      assert( discriminant >= 0 );
      NT d = sqrt(discriminant);
      if ((smaller && a>0) || (!smaller && a<0))
        d = -d;
      return (d-b)/(a*2);
    }

    // This version is internal and can be re-used for
    // number types which are rational.
    template < typename RT, typename FT >
    Root_of_2< RT >
    make_root_of_2_rational(const FT &a, const FT &b, const FT &c, bool smaller)
    {
      typedef CGAL::Rational_traits< FT > Rational;

      Rational r;
      assert( r.denominator(a) > 0 );
      assert( r.denominator(b) > 0 );
      assert( r.denominator(c) > 0 );

/*   const RT lcm = ( r.denominator(a) * r.denominator(b) * r.denominator(c)          )/
               ( gcd( r.denominator(a), gcd(r.denominator(b), r.denominator(c)) ) );

      RT a_ = r.numerator(a) * ( lcm / r.denominator(a) );
      RT b_ = r.numerator(b) * ( lcm / r.denominator(b) );
      RT c_ = r.numerator(c) * ( lcm / r.denominator(c) );
*/
      RT a_ = r.numerator(a) * r.denominator(b) * r.denominator(c);
      RT b_ = r.numerator(b) * r.denominator(a) * r.denominator(c);
      RT c_ = r.numerator(c) * r.denominator(a) * r.denominator(b);

      return make_root_of_2(a_,b_,c_,smaller);
    }

  } // namespace CGALi

template < typename RT >
Sign
sign(const Root_of_2<RT> &a)
{
  // We use an optimized version of the equivalent to :
  // return static_cast<Sign>((int) compare(a, 0));

  assert(is_valid(a));
  // First, we compare 0 to the root of the derivative of a.
  // (which is equivalent to a[1])

  int sgn = sign(a[1]);

  if (sgn > 0)
    return a.is_smaller() ? NEGATIVE : opposite(sign(a[0]));

  if (sgn < 0)
    return a.is_smaller() ? sign(a[0]) : POSITIVE;

  if (is_zero(a[0]))
    return ZERO;
  return a.is_smaller() ? NEGATIVE : POSITIVE;
}


template < typename RT >
Comparison_result
compare(const Root_of_2<RT> &a, const Root_of_2<RT> &b)
{
  std::cout<< "a: " << CGAL::to_double(a) 
	   << " b: " << CGAL::to_double(b) << std::endl;
  assert(is_valid(a) && is_valid(b));

  // Now a and b are both of degree 2.
  if (a.is_smaller())
    {
      if (b.is_smaller())
          return CGALi::compare_22_11(a[2], a[1], a[0], b[2], b[1], b[0]);

      return CGALi::compare_22_12(a[2], a[1], a[0], b[2], b[1], b[0]);
    }
  if (b.is_smaller())
    return CGALi::compare_22_21(a[2], a[1], a[0], b[2], b[1], b[0]);

  return CGALi::compare_22_22(a[2], a[1], a[0], b[2], b[1], b[0]);
}

template < typename RT >
Comparison_result
compare(const Root_of_2<RT> &a,
	const typename Root_of_traits< RT >::RootOf_1 &b)
{
  typedef typename Root_of_traits< RT >::RootOf_1 RootOf_1;
  
  assert(is_valid(a) && is_valid(b));

  RootOf_1 d_a(-a[1],2*a[2]);

  int cmp = compare(b,d_a);

  if (cmp > 0)
    return (a.is_smaller() ?
      SMALLER :
      static_cast<Comparison_result>( -a.sign_at(b)));

  if (cmp < 0)
    return (a.is_smaller() ?
      static_cast<Comparison_result>((int) a.sign_at(b)) :
      LARGER);

  if (is_zero(a.discriminant()))
    return EQUAL;
  return (a.is_smaller() ? SMALLER : LARGER);
}

template < typename RT >  inline
Comparison_result
compare(const typename Root_of_traits< RT >::RootOf_1 &a,
	const Root_of_2<RT> &b)
{
   return opposite(compare(b, a));
}

template < typename RT >
Comparison_result
compare(const Root_of_2<RT> &a, const RT &b)
{
  assert(is_valid(a));

  // First, we compare b to the root of the derivative of a.
  int cmp = compare(2*a[2]*b, -a[1]);

  if (cmp > 0)
    return a.is_smaller() ? SMALLER
                          : (Comparison_result) - a.sign_at(b);
  if (cmp < 0)
    return a.is_smaller() ? (Comparison_result) a.sign_at(b)
                          : LARGER;

  if (is_zero(a.discriminant()))
    return EQUAL;
  return a.is_smaller() ? SMALLER : LARGER;
}

template < typename RT >  inline
Comparison_result
compare(const RT &a, const Root_of_2<RT> &b)
{
   return opposite(compare(b, a));
}

template < typename RT >  inline
Comparison_result
compare(const Root_of_2<RT> &a, const CGAL_CK_int(RT) &b)
{
   return compare(a, RT(b));
}

template < typename RT >  inline
Comparison_result
compare(const CGAL_CK_int(RT) &a, const Root_of_2<RT> &b)
{
   return opposite(compare(b, RT(a)));
}


template < typename RT >
inline
bool
operator<(const Root_of_2<RT> &a, const Root_of_2<RT> &b)
{
  return compare(a, b) < 0;
}

template < typename RT >
inline
bool
operator<(const typename Root_of_traits< RT >::RootOf_1 &a,
	  const Root_of_2<RT> &b)
{
  return compare(a, b) < 0;
}

template < typename RT >
inline
bool
operator<(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1 &b)
{
  return compare(a, b) < 0;
}

template < typename RT >
inline
bool
operator<(const RT &a, const Root_of_2<RT> &b)
{
  return compare(a, b) < 0;
}

template < typename RT >
inline
bool
operator<(const Root_of_2<RT> &a, const RT &b)
{
  return compare(a, b) < 0;
}

template < typename RT >
inline
bool
operator<(const CGAL_CK_int(RT) &a, const Root_of_2<RT> &b)
{
  return compare(a, b) < 0;
}

template < typename RT >
inline
bool
operator<(const Root_of_2<RT> &a, const CGAL_CK_int(RT) &b)
{
  return compare(a, b) < 0;
}


template < typename RT >
inline
bool
operator>(const Root_of_2<RT> &a, const Root_of_2<RT> &b)
{
  return b < a;
}

template < typename RT >
inline
bool
operator>(const typename Root_of_traits< RT >::RootOf_1 &a,
	  const Root_of_2<RT> &b)
{
  return b < a;
}

template < typename RT >
inline
bool
operator>(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1 &b)
{
  return b < a;
}

template < typename RT >
inline
bool
operator>(const RT &a, const Root_of_2<RT> &b)
{
  return b < a;
}

template < typename RT >
inline
bool
operator>(const Root_of_2<RT> &a, const RT &b)
{
  return b < a;
}

template < typename RT >
inline
bool
operator>(const CGAL_CK_int(RT) &a, const Root_of_2<RT> &b)
{
  return b < a;
}

template < typename RT >
inline
bool
operator>(const Root_of_2<RT> &a, const CGAL_CK_int(RT) &b)
{
  return b < a;
}


template < typename RT >
inline
bool
operator>=(const Root_of_2<RT> &a, const Root_of_2<RT> &b)
{
  return !(a < b);
}

template < typename RT >
inline
bool
operator>=(const typename Root_of_traits< RT >::RootOf_1 &a,
	   const Root_of_2<RT> &b)
{
  return !(a < b);
}

template < typename RT >
inline
bool
operator>=(const Root_of_2<RT> &a,
	   const typename Root_of_traits< RT >::RootOf_1 &b)
{
  return !(a < b);
}

template < typename RT >
inline
bool
operator>=(const RT &a, const Root_of_2<RT> &b)
{
  return !(a < b);
}

template < typename RT >
inline
bool
operator>=(const Root_of_2<RT> &a, const RT &b)
{
  return !(a < b);
}

template < typename RT >
inline
bool
operator>=(const CGAL_CK_int(RT) &a, const Root_of_2<RT> &b)
{
  return !(a < b);
}

template < typename RT >
inline
bool
operator>=(const Root_of_2<RT> &a, const CGAL_CK_int(RT) &b)
{
  return !(a < b);
}


template < typename RT >
inline
bool
operator<=(const Root_of_2<RT> &a, const Root_of_2<RT> &b)
{
  return !(a > b);
}

template < typename RT >
inline
bool
operator<=(const typename Root_of_traits< RT >::RootOf_1 &a,
	   const Root_of_2<RT> &b)
{
  return !(a > b);
}

template < typename RT >
inline
bool
operator<=(const Root_of_2<RT> &a,
	   const typename Root_of_traits< RT >::RootOf_1 &b)
{
  return !(a > b);
}

template < typename RT >
inline
bool
operator<=(const RT &a, const Root_of_2<RT> &b)
{
  return !(a > b);
}

template < typename RT >
inline
bool
operator<=(const Root_of_2<RT> &a, const RT &b)
{
  return !(a > b);
}

template < typename RT >
inline
bool
operator<=(const CGAL_CK_int(RT) &a, const Root_of_2<RT> &b)
{
  return !(a > b);
}

template < typename RT >
inline
bool
operator<=(const Root_of_2<RT> &a, const CGAL_CK_int(RT) &b)
{
  return !(a > b);
}


template < typename RT >
inline
bool
operator==(const Root_of_2<RT> &a, const Root_of_2<RT> &b)
{
  return compare(a, b) == 0;
}

template < typename RT >
inline
bool
operator==(const typename Root_of_traits< RT >::RootOf_1 &a,
	   const Root_of_2<RT> &b)
{
  return compare(a, b) == 0;
}

template < typename RT >
inline
bool
operator==(const Root_of_2<RT> &a,
	   const typename Root_of_traits< RT >::RootOf_1 &b)
{
  return compare(a, b) == 0;
}

template < typename RT >
inline
bool
operator==(const RT &a, const Root_of_2<RT> &b)
{
  return compare(a, b) == 0;
}

template < typename RT >
inline
bool
operator==(const Root_of_2<RT> &a, const RT &b)
{
  return compare(a, b) == 0;
}

template < typename RT >
inline
bool
operator==(const CGAL_CK_int(RT) &a, const Root_of_2<RT> &b)
{
  return compare(a, b) == 0;
}

template < typename RT >
inline
bool
operator==(const Root_of_2<RT> &a, const CGAL_CK_int(RT) &b)
{
  return compare(a, b) == 0;
}


template < typename RT >
inline
bool
operator!=(const Root_of_2<RT> &a, const Root_of_2<RT> &b)
{
  return !(a == b);
}

template < typename RT >
inline
bool
operator!=(const typename Root_of_traits< RT >::RootOf_1 &a,
	   const Root_of_2<RT> &b)
{
  return !(a == b);
}

template < typename RT >
inline
bool
operator!=(const Root_of_2<RT> &a,
	   const typename Root_of_traits< RT >::RootOf_1 &b)
{
  return !(a == b);
}

template < typename RT >
inline
bool
operator!=(const RT &a, const Root_of_2<RT> &b)
{
  return !(a == b);
}

template < typename RT >
inline
bool
operator!=(const Root_of_2<RT> &a, const RT &b)
{
  return !(a == b);
}

template < typename RT >
inline
bool
operator!=(const CGAL_CK_int(RT) &a, const Root_of_2<RT> &b)
{
  return !(a == b);
}

template < typename RT >
inline
bool
operator!=(const Root_of_2<RT> &a, const CGAL_CK_int(RT) &b)
{
  return !(a == b);
}


template < typename RT >
Root_of_2<RT>
square(const Root_of_2<RT> &a)
{
  assert(is_valid(a));

  // It's easy to get the explicit formulas for the square of the two roots.
  // Then it's easy to compute their sum and their product, which gives the
  // coefficients of the polynomial (X^2 - Sum X + Product).
 return Root_of_2<RT> ( square(a[2]),
                        2 * a[2] * a[0] - square(a[1]),
                        square(a[0]),
                        (a.is_smaller() ^ (a[1]>0)));
}

template < typename RT >
Root_of_2<RT>
operator-(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1& b)
{
  typedef typename Root_of_traits< RT >::RootOf_1  RO1;
  typedef CGAL::Rational_traits< RO1 >             Rational;
  //RT should be the same as Rational::RT

  Rational r;
  assert(is_valid(a) && is_valid(b));
  const RT &b1 = r.denominator(b);
  const RT &b0 = r.numerator(b);
  //assert(b1>0);

  RT sqb1 = square(b1);
  RT sqb0 = square(b0);
  RT b0b1 = b0 * b1;

  return Root_of_2<RT>(a[2] * sqb1,
                       2*a[2]*b0b1 + a[1]*sqb1,
                       a[2]* sqb0 + a[1]*b0b1 + a[0]*sqb1,
                       a.is_smaller());
}

template < typename RT >
inline
Root_of_2<RT>
operator-(const typename Root_of_traits< RT >::RootOf_1 &a,
	  const Root_of_2<RT> &b)
{
  return -(b-a);
}

template < typename RT >
Root_of_2<RT>
operator-(const Root_of_2<RT> &a, const RT& b)
{
  assert(is_valid(a) && is_valid(b));

  // It's easy to see it using the formula (X^2 - sum X + prod).

  //RT p = a[0] + b * a[1] + a[2] * square(b);
  //RT s = a[1] + b * a[2] * 2;

  RT tmp = a[2] * b;
  RT s = a[1] + tmp;
  RT p = a[0] + s * b;
  s += tmp;

  return Root_of_2<RT>(a[2], s, p, a.is_smaller());
}

template < typename RT >
inline
Root_of_2<RT>
operator-(const RT &a, const Root_of_2<RT> &b)
{
  return -(b-a);
}

template < typename RT >
inline
Root_of_2<RT>
operator+(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1& b)
{
  return a - typename Root_of_traits< RT >::RootOf_1(-b);
}

template < typename RT >
inline
Root_of_2<RT>
operator+(const typename Root_of_traits< RT >::RootOf_1 &a,
	  const Root_of_2<RT> &b)
{
  return b - typename Root_of_traits< RT >::RootOf_1(-a);
}

template < typename RT >
inline
Root_of_2<RT>
operator+(const Root_of_2<RT> &a, const RT& b)
{
  return a - RT(-b);
}

template < typename RT >
inline
Root_of_2<RT>
operator+(const RT &a, const Root_of_2<RT> &b)
{
  return b - RT(-a);
}

template < typename RT >
Root_of_2<RT>
operator*(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1& b)
{
  typedef typename Root_of_traits< RT >::RootOf_1  RO1;
  typedef CGAL::Rational_traits< RO1 >             Rational;
  //RT should be the same as Rational::RT

  Rational r;
  assert(is_valid(a) && is_valid(b));
  const RT &b1 = r.denominator(b);
  const RT &b0 = r.numerator(b);
  assert(b1>0);

  return Root_of_2<RT>(a[2] * square(b1), a[1] * b1 * b0,
                       a[0] * square(b0),
                       sign(b) < 0 ? !a.is_smaller()
                                   :  a.is_smaller());
}

template < typename RT >
inline
Root_of_2<RT>
operator*(const typename Root_of_traits< RT >::RootOf_1 &a,
	  const Root_of_2<RT> &b)
{
  return b * a;
}

template < typename RT >
Root_of_2<RT>
operator*(const Root_of_2<RT> &a, const RT& b)
{
  assert(is_valid(a));

  return Root_of_2<RT>(a[2], a[1] * b, a[0] * square(b),
                         b < 0 ? !a.is_smaller() : a.is_smaller());
}

template < typename RT >
inline
Root_of_2<RT>
operator*(const RT &a, const Root_of_2<RT> &b)
{
  return b * a;
}


template < typename RT >
double
to_double(const Root_of_2<RT> &x)
{
  assert(is_valid(x));

  double a = to_double(x[2]);
  double b = to_double(x[1]);
  double d = std::sqrt(to_double(x.discriminant()));

  assert(a > 0);
  if (x.is_smaller())
    d = -d;

  return (d-b)/(a*2);
}

template < typename RT >
std::pair<double, double>
to_interval(const Root_of_2<RT> &x)
{
  assert(is_valid(x));

  CGAL::Interval_nt<> a = to_interval(x[2]);
  CGAL::Interval_nt<> b = to_interval(x[1]);
  CGAL::Interval_nt<> disc = to_interval(x.discriminant());
  CGAL::Interval_nt<> d = sqrt(disc);

  if (x.is_smaller())
    d = -d;

  return ((d-b)/(a*2)).pair();
}

template < typename RT >
bool
is_valid(const Root_of_2<RT> &r)
{
  return r.is_valid();
}

template < typename RT >
std::ostream &
operator<<(std::ostream &os, const Root_of_2<RT> &r)
{
  return os << to_double(r);
}

template < typename RT >
void
print(std::ostream &os, const Root_of_2<RT> &r)
{
    os << "Root_of_2( (" 
       << r[2] << ") * X^2 + ("
       << r[1] << ") * X + ("
       << r[0] << ") , "
       << (r.is_smaller() ? "smaller" : "larger") << " )";
}


} // namespace CGAL


// CGAL requirements for number types.

#include <CGAL/tags.h>

namespace CGAL {

using CGAL::is_valid;
using CGAL::to_double;
using CGAL::to_interval;
using CGAL::compare;
using CGAL::sign;
using CGAL::square;

template < class RT >
struct Number_type_traits < CGAL::Root_of_2<RT> > {
  typedef Tag_false      Has_gcd;
  typedef Tag_false      Has_division;
  typedef Tag_false      Has_sqrt;
};

template < typename RT >
inline
io_Operator
io_tag(const CGAL::Root_of_2<RT>&)
{
  return io_Operator();
}

} // namespace CGAL

#endif // CGAL_ROOT_OF_ROOT_OF_2_H
