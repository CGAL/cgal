// Copyright (c) 2005,2006  INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion, Monique Teillaud, Athanasios Kakargias

#ifndef CGAL_ROOT_OF_2_H
#define CGAL_ROOT_OF_2_H

#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Root_of_2_fwd.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/NT_converter.h>

#include <CGAL/enum.h>
#include <CGAL/tags.h>
#include <CGAL/Number_types/internal_functions_comparison_root_of_2.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Quotient.h>
#include <CGAL/assertions.h>
#include <CGAL/Binary_operator_result.h>
#include <boost/type_traits/is_same.hpp>

#define CGAL_int(T)    typename First_if_different<int,    T>::Type
#define CGAL_double(T) typename First_if_different<double, T>::Type

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
// - add sqrt() (when it's degree 1), or a make_sqrt<RT>(const RT &r) ?
// - inverse()
// TODO :
// - use Boost.Operators.
// - add subtraction/addition with a degree 2 Root_of of the same field ?
// - add +, -, *, / (when it is possible) ?
// - add constructor from CGAL::Polynomial ?
//   There should be a proper separate class Polynomial.
// - in compare_roots, we evaluate the polynomial at some FT, or at some
//   root of degree 1 polynomials.  It would be nice to have a separate
//   polynomial class which performed this task (and others)...
// - overloaded versions of make_root_of_2<>() for Lazy_exact_nt<> and others.

template < typename RT_ >
class Root_of_2 {

  RT_  C0,C1,C2;   // Coefficients (see below)
                           // _smaller indicates if it's the smaller of the 2 roots.
  unsigned char _smaller;  // (we use a unsigned char because it's smaller than a bool)
  unsigned char _rational;
  unsigned char _delta_is_not_zero;

  // the value is the root of P(X) = C2.X^2 + C1.X + C0,
  // and C2 > 0.

public:

  typedef RT_    RT;
  typedef typename Root_of_traits<RT>::RootOf_1  FT;

  Root_of_2()
    : C0(0), C1(0), C2(1), _rational(1)
  {
    CGAL_assertion(is_valid());
  }

  Root_of_2(const RT& c0)
    : C1(c0), C2(1), _rational(1)
  {
    CGAL_assertion(is_valid());
  }

  Root_of_2(const typename First_if_different<int, RT>::Type & c0)
    : C1(RT(c0)), C2(1), _rational(1)
  {
    CGAL_assertion(is_valid());
  }

  Root_of_2(const typename First_if_different<FT, RT>::Type & c0)
  {
    typedef CGAL::Rational_traits< FT > Rational;

    Rational r;
    CGAL_assertion(!is_zero(r.denominator(c0)));
    *this = Root_of_2(r.numerator(c0), r.denominator(c0));

    CGAL_assertion(is_valid());
  }

  Root_of_2(const RT& c1, const RT& c0)
    : C1(c1), C2(c0), _rational(1)
  {
    // It reduces the performance and we don't seem to need it
    // simplify_quotient(C1, C2);  
    CGAL_assertion(is_valid());
  }

  Root_of_2(const RT& a, const RT& b, const RT& c, const bool s, 
            const bool dinz = false)
    : C0(c), C1(b), C2(a), _smaller(s), _rational(0), 
      _delta_is_not_zero(dinz)
  {
    CGAL_assertion(a != 0);
    if (a < 0) {
      C0 = -c;
      C1 = -b;
      C2 = -a;
    }
    simplify_root_of_2(C2,C1,C0);
    CGAL_assertion(is_valid());
  }
  
  /* Jetter un coup d'oeil pour voir s'il n'y a aucun moyen de faire comme ca

  template <typename RT2>
  Root_of_2(const Root_of_2<RT2>& r)
  {
    CGAL_assertion(is_valid());
    if(r.is_rational()) {
      C1 = r[1]; 
      C2 = r[2];
      _rational = 1;
    } else {
      C0 = r[0];
      C1 = r[1]; 
      C2 = r[2];
      _smaller = r.is_smaller();
      _rational = 0;
      _delta_is_not_zero = r.is_known_delta_not_zero();
    }
  }*/

  template <typename RT2>
  Root_of_2(const Root_of_2<RT2>& r)
    : C0(r[0]), C1(r[1]), C2(r[2]), _smaller(r.is_smaller()),
      _delta_is_not_zero(r.is_known_delta_not_zero())
  {
    CGAL_assertion(is_valid());
  }

  Root_of_2 operator-() const
  {
    if(is_rational()) return Root_of_2(-C1,C2);
    return Root_of_2 (C2, -C1, C0, !is_smaller(), 
      is_known_delta_not_zero());
  }

  bool is_valid() const
  {
    if(is_rational()) {
      return !is_zero(C2);
    } else {
      return (C2 > 0) && ! is_negative(discriminant());
    }
  }

  bool is_smaller() const
  {
    return _smaller;
  }

  bool is_rational() const
  {
    return _rational;
  }

  bool is_known_delta_not_zero() const 
  {
    return _delta_is_not_zero;
  }

  Root_of_2 inverse() const 
  {
    if(is_rational()) {
      CGAL_assertion(C1 != 0);
      return Root_of_2(C2,C1);
    }
    if(is_zero(C0)) {
	CGAL_assertion(C1 != 0);
        CGAL_assertion((C1 > 0 && _smaller) || (C1 < 0 && !_smaller));
        return Root_of_2(-C2,C1);
    } else if(C0 > 0) return Root_of_2(C0, C1, C2, !_smaller,
      is_known_delta_not_zero());
    else return Root_of_2(C0, C1, C2, _smaller,
      is_known_delta_not_zero());
  }

  // The following functions deal with the internal polynomial.
  // Probably they should move to a separate polynomial class.

  const RT & operator[](int i) const
  {
    CGAL_assertion(i<3);
    return (&C0)[i];
  }

  RT discriminant() const
  {
    if(is_rational()) return RT(0);
    return CGAL_NTS square(C1) - 4*C0*C2;
  }

  template < typename T > 
  T eval_at(const T& x) const 
  {
    if(is_rational()) return x * C2 - C1;
    if(is_zero(x)) return C0;
    const bool zeroC0 = is_zero(C0);
    const bool zeroC1 = is_zero(C1);
    if(zeroC0 && zeroC1) return x * C2;
    if(zeroC0) return x * (C1 + x * C2);
    if(zeroC1) return (x * x * C2) + C0;
    return C0 + x * (C1 + x * C2); 
  } 

  template < typename T > 
  Sign sign_at(const T &x) const
  {
    // Maybe there is slightly more efficient.
    return CGAL_NTS sign(eval_at(x));
  }

  Root_of_2 conjugate() const
  {
    if(is_rational()) return Root_of_2(C1, C2);
    return Root_of_2(C2, C1, C0, !is_smaller(), 
      is_known_delta_not_zero());
  }
}; // Root_of_2


template < class NT1,class NT2 >
struct NT_converter < Root_of_2<NT1> , Root_of_2<NT2> >
  : public std::unary_function< NT1, NT2 >
{
    Root_of_2<NT2>
    operator()(const Root_of_2<NT1> &a) const
    {
      if(!a.is_rational()) {
        return make_root_of_2(NT_converter<NT1,NT2>()(a[2]),NT_converter<NT1,NT2>()(a[1]),
                              NT_converter<NT1,NT2>()(a[0]),a.is_smaller());
      } else {
        return make_root_of_2(NT_converter<NT1,NT2>()(a[1]) / NT_converter<NT1,NT2>()(a[2]));
      }
    }
};

template < class NT1,class NT2 >
struct NT_converter < NT1 , Root_of_2<NT2> >
  : public std::unary_function< NT1, NT2 >
{
    Root_of_2<NT2>
    operator()(const NT1 &a) const
    {
        return Root_of_2<NT2>(NT_converter<NT1,NT2>()(a));
    }
};





template < class NT1 >
struct NT_converter < Root_of_2<NT1>, Root_of_2<NT1> >
  : public std::unary_function< NT1, NT1 >
{
    const Root_of_2<NT1> &
    operator()(const Root_of_2<NT1> &a) const
    {
        return a;
    }
};

// Simplify the coefficients of root_of_2.
// Currently the default template doesn't do anything.
// This function is not documented as a number type requirement for now.
template < typename RT >
inline void
simplify_root_of_2(RT &, RT &, RT&) {}

template < typename RT >
Sign
sign(const Root_of_2<RT> &a)
{

  CGAL_assertion(is_valid(a));

  if(a.is_rational()) {
    if(is_zero(a[1])) return ZERO;
    if(sign(a[1]) == sign(a[2])) return POSITIVE;
    return NEGATIVE;
  }

  // We use an optimized version of the equivalent to :
  // return static_cast<Sign>((int) compare(a, 0));

  // First, we compare 0 to the root of the derivative of a.
  // (which is equivalent to a[1])

  int sgn = CGAL_NTS sign(a[1]);

  if (sgn > 0)
    return a.is_smaller() ? NEGATIVE : opposite(CGAL_NTS sign(a[0]));

  if (sgn < 0)
    return a.is_smaller() ? CGAL_NTS sign(a[0]) : POSITIVE;

  if (CGAL_NTS is_zero(a[0]))
    return ZERO;
  return a.is_smaller() ? NEGATIVE : POSITIVE;
}

// A namespace for some internal comparison functions.
// The idea is that those functions are more optimized then
// the generics one, and it is used only for Root_of_2 comparisons
namespace CGALi {

// Compare xnum/xden with n
template <class NT>
inline
Comparison_result
compare(const NT & xnum, const NT & xden,
             const NT & n) 
{
  Sign sig_xnum = CGAL_NTS sign(xnum);
  Sign sig_xden = CGAL_NTS sign(xden);
  Sign ysign = CGAL_NTS sign(n);
  int xsign = sig_xnum * sig_xden;
  if (xsign == 0) return static_cast<Comparison_result>(-ysign);
  if (ysign == 0) return static_cast<Comparison_result>(xsign);
  int diff = xsign - ysign;
  if (diff == 0) {
    if(sig_xden > 0) return CGAL_NTS compare(xnum, n * xden);
    return opposite(CGAL_NTS compare(xnum, n * xden)); 
  } else {
    return (xsign < ysign) ? SMALLER : LARGER;
  }
}

// compare xnum/xden with ynum/yden
template <class NT>
inline
Comparison_result
compare(const NT & xnum, const NT & xden,
             const NT & ynum, const NT & yden) 
{
  Sign sig_xnum = CGAL_NTS sign(xnum);
  Sign sig_xden = CGAL_NTS sign(xden);
  Sign sig_ynum = CGAL_NTS sign(ynum);
  Sign sig_yden = CGAL_NTS sign(yden);
  int xsign = sig_xnum * sig_xden;
  int ysign = sig_ynum * sig_yden;
  if (xsign == 0) return static_cast<Comparison_result>(-ysign);
  if (ysign == 0) return static_cast<Comparison_result>(xsign);
  int diff = xsign - ysign;
  if (diff == 0) {
    if((sig_xden * sig_yden) > 0) 
      return CGAL_NTS compare(xnum * yden, ynum * xden);
    return opposite(CGAL_NTS compare(xnum * yden, ynum * xden));
  } else {
    return (xsign < ysign) ? SMALLER : LARGER;
  }
}

// with a = Ax^2 + Bx + C
// sign(a.evaluate(x)) with x = FT (n1/n2).
// equivalent to opposite(compare(An1^2 + Cn2^2, Bn1n2))
template < typename RT >
inline
Comparison_result
evaluate_at(const Root_of_2<RT> &a, const RT & bnum, const RT & bden) {
  return CGAL_NTS compare( 
    bnum*(a[2] * bnum) + bden * (a[0] * bden),
    - (a[1] * bnum) * bden 
  );
}

// with a = Ax^2 + Bx + C
// -sign(a.evaluate(x)) with x = FT (n1/n2).
// equivalent to compare(An1^2 + Cn2^2, Bn1n2)
template < typename RT >
inline
Comparison_result
evaluate_neg_at(const Root_of_2<RT> &a, const RT & bnum, const RT & bden) {
  return opposite(evaluate_at(a,bnum,bden));
}

// compare Root_of_2<RT> a with bnum/bden
template < typename RT >
inline
Comparison_result
compare(const Root_of_2<RT> &a, const RT & bnum, const RT & bden) 
{
  typedef typename Root_of_traits< RT >::RootOf_1 RootOf_1;

  int cmp = CGALi::compare(bnum,bden,-a[1],2*a[2]);

  if(is_zero(bnum)) {
    if (cmp > 0) return (a.is_smaller() ? SMALLER : 
      static_cast<Comparison_result>(-sign(a[0])));
    if (cmp < 0) return (a.is_smaller() ? 
      static_cast<Comparison_result>((int) sign(a[0])) : LARGER);
  } 

  if (cmp > 0)
    return (a.is_smaller() ?
      SMALLER :
      CGALi::evaluate_neg_at(a,bnum,bden));

  if (cmp < 0)
    return (a.is_smaller() ?
      CGALi::evaluate_at(a,bnum,bden) :
      LARGER);

  // This case is rational, whenever we have 
  // a rational case we should use the rational case constructor
  // whenever we are sure that the delta is different from 0, put on constructor
  if(a.is_known_delta_not_zero()) 
    return (a.is_smaller() ? SMALLER : LARGER);
  if (is_zero(a.discriminant())) return EQUAL;
  return (a.is_smaller() ? SMALLER : LARGER);
}

}

template < typename RT >
Comparison_result
compare(const Root_of_2<RT> &a, const Root_of_2<RT> &b)
{
  typedef typename Root_of_traits< RT >::RootOf_1 RootOf_1;
  typedef typename First_if_different<RootOf_1, RT>::Type WhatType;
  typedef typename boost::is_same< WhatType, RT > do_not_filter;

  CGAL_assertion(is_valid(a) && is_valid(b));

  if(!do_not_filter::value) {
    Interval_nt<> ia = to_interval(a);
    Interval_nt<> ib = to_interval(b);
    if(ia.inf() > ib.sup()) return LARGER;
    if(ia.sup() < ib.inf()) return SMALLER;
  }

  if(a.is_rational() && b.is_rational()) {
    return CGALi::compare(a[1], a[2], b[1], b[2]);
  }

  if(a.is_rational()) {
    return opposite(CGAL_NTS CGALi::compare(b, a[1], a[2]));
  }

  if(b.is_rational()) {
    return CGALi::compare(a, b[1], b[2]);
  }

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
  typedef typename First_if_different<RootOf_1, RT>::Type WhatType;
  typedef typename boost::is_same< WhatType, RT > do_not_filter;

  CGAL_assertion(is_valid(a) && is_valid(b));

  if(!do_not_filter::value) {
    Interval_nt<> ia = to_interval(a);
    Interval_nt<> ib = to_interval(b);
    if(ia.inf() > ib.sup()) return LARGER;
    if(ia.sup() < ib.inf()) return SMALLER;
  }

  if(a.is_rational()) {
    CGAL::Rational_traits< RootOf_1 > r;
    return CGALi::compare(a[1], a[2], r.numerator(b), r.denominator(b));
  }

  RootOf_1 d_a(-a[1],2*a[2]);

  int cmp = CGAL_NTS compare(b,d_a);

  if(is_zero(b)) {
    if (cmp > 0) return (a.is_smaller() ? SMALLER : 
      static_cast<Comparison_result>(-sign(a[0])));
    if (cmp < 0) return (a.is_smaller() ? 
      static_cast<Comparison_result>((int) sign(a[0])) : LARGER);
  } 

  CGAL::Rational_traits< RootOf_1 > r;

  if (cmp > 0)
    return (a.is_smaller() ?
      SMALLER :
      CGALi::evaluate_neg_at(a, r.numerator(b), r.denominator(b)));

  if (cmp < 0)
    return (a.is_smaller() ?
      CGALi::evaluate_at(a, r.numerator(b), r.denominator(b)) :
      LARGER);

  if(a.is_known_delta_not_zero()) 
    return (a.is_smaller() ? SMALLER : LARGER);
  if (is_zero(a.discriminant())) return EQUAL;
  return (a.is_smaller() ? SMALLER : LARGER);
}

template < typename RT >  inline
Comparison_result
compare(const typename Root_of_traits< RT >::RootOf_1 &a,
	const Root_of_2<RT> &b)
{
   return opposite(CGAL_NTS compare(b, a));
}

template < typename RT >
Comparison_result
compare(const Root_of_2<RT> &a, const RT &b)
{
  typedef typename Root_of_traits< RT >::RootOf_1 RootOf_1;
  typedef typename First_if_different<RootOf_1, RT>::Type WhatType;
  typedef typename boost::is_same< WhatType, RT > do_not_filter;

  CGAL_assertion(is_valid(a) && is_valid(b));

  if(!do_not_filter::value) {
    Interval_nt<> ia = to_interval(a);
    Interval_nt<> ib = to_interval(b);
    if(ia.inf() > ib.sup()) return LARGER;
    if(ia.sup() < ib.inf()) return SMALLER;
  }

  if(a.is_rational()) {
    return CGALi::compare(a[1], a[2], b);
  }

  // First, we compare b to the root of the derivative of a.
  int cmp = CGAL_NTS compare(2*a[2]*b, -a[1]);

  if(is_zero(b)) {
    if (cmp > 0) return (a.is_smaller() ? SMALLER : 
      static_cast<Comparison_result>(-sign(a[0])));
    if (cmp < 0) return (a.is_smaller() ? 
      static_cast<Comparison_result>((int) sign(a[0])) : LARGER);
  } 

  if (cmp > 0)
    return (a.is_smaller() ?
      SMALLER :
      static_cast<Comparison_result>( -a.sign_at(b)));

  if (cmp < 0)
    return (a.is_smaller() ?
      static_cast<Comparison_result>((int) a.sign_at(b)) :
      LARGER);

  if(a.is_known_delta_not_zero()) 
    return (a.is_smaller() ? SMALLER : LARGER);
  if (is_zero(a.discriminant())) return EQUAL;
  return (a.is_smaller() ? SMALLER : LARGER);
}

template < typename RT >  inline
Comparison_result
compare(const RT &a, const Root_of_2<RT> &b)
{
   return opposite(CGAL_NTS compare(b, a));
}

template < typename RT >  inline
Comparison_result
compare(const Root_of_2<RT> &a, const CGAL_int(RT) &b)
{
   return CGAL_NTS compare(a, RT(b));
}

template < typename RT >  inline
Comparison_result
compare(const CGAL_int(RT) &a, const Root_of_2<RT> &b)
{
   return opposite(CGAL_NTS compare(b, RT(a)));
}


template < typename RT > inline
bool operator<(const Root_of_2<RT> &a, const Root_of_2<RT> &b) {
  return CGAL_NTS compare(a, b) < 0;
}

template < typename RT > inline
bool operator<(const typename Root_of_traits< RT >::RootOf_1 &a,
	  const Root_of_2<RT> &b) 
{
  return CGAL_NTS compare(a, b) < 0;
}

template < typename RT > inline
bool operator<(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1 &b) 
{
  return CGAL_NTS compare(a, b) < 0;
}

template < typename RT > inline
bool operator<(const RT &a, const Root_of_2<RT> &b) 
{
  return CGAL_NTS compare(a, b) < 0;
}

template < typename RT > inline
bool operator<(const Root_of_2<RT> &a, const RT &b) 
{
  return CGAL_NTS compare(a, b) < 0;
}

template < typename RT > inline
bool operator<(const CGAL_int(RT) &a, const Root_of_2<RT> &b) 
{
  return CGAL_NTS compare(a, b) < 0;
}

template < typename RT > inline
bool operator<(const Root_of_2<RT> &a, const CGAL_int(RT) &b) 
{
  return CGAL_NTS compare(a, b) < 0;
}

template < typename RT > inline
bool operator>(const Root_of_2<RT> &a, const Root_of_2<RT> &b) 
{
  return b < a;
}

template < typename RT > inline
bool operator>(const typename Root_of_traits< RT >::RootOf_1 &a,
	  const Root_of_2<RT> &b) 
{
  return b < a;
}

template < typename RT > inline
bool operator>(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1 &b) 
{
  return b < a;
}

template < typename RT > inline
bool operator>(const RT &a, const Root_of_2<RT> &b) 
{
  return b < a;
}

template < typename RT > inline
bool operator>(const Root_of_2<RT> &a, const RT &b) 
{
  return b < a;
}

template < typename RT > inline
bool operator>(const CGAL_int(RT) &a, const Root_of_2<RT> &b) 
{
  return b < a;
}

template < typename RT > inline
bool operator>(const Root_of_2<RT> &a, const CGAL_int(RT) &b) 
{
  return b < a;
}

template < typename RT > inline
bool operator>=(const Root_of_2<RT> &a, const Root_of_2<RT> &b) 
{
  return !(a < b);
}

template < typename RT > inline
bool operator>=(const typename Root_of_traits< RT >::RootOf_1 &a,
	   const Root_of_2<RT> &b) 
{
  return !(a < b);
}

template < typename RT > inline
bool operator>=(const Root_of_2<RT> &a,
	   const typename Root_of_traits< RT >::RootOf_1 &b) 
{
  return !(a < b);
}

template < typename RT > inline
bool operator>=(const RT &a, const Root_of_2<RT> &b) 
{
  return !(a < b);
}

template < typename RT > inline
bool operator>=(const Root_of_2<RT> &a, const RT &b) 
{
  return !(a < b);
}

template < typename RT > inline
bool operator>=(const CGAL_int(RT) &a, const Root_of_2<RT> &b) 
{
  return !(a < b);
}

template < typename RT > inline
bool operator>=(const Root_of_2<RT> &a, const CGAL_int(RT) &b) 
{
  return !(a < b);
}

template < typename RT > inline
bool operator<=(const Root_of_2<RT> &a, const Root_of_2<RT> &b) 
{
  return !(a > b);
}

template < typename RT > inline
bool operator<=(const typename Root_of_traits< RT >::RootOf_1 &a,
	   const Root_of_2<RT> &b) 
{
  return !(a > b);
}

template < typename RT > inline
bool operator<=(const Root_of_2<RT> &a,
	   const typename Root_of_traits< RT >::RootOf_1 &b) 
{
  return !(a > b);
}

template < typename RT > inline
bool operator<=(const RT &a, const Root_of_2<RT> &b) 
{
  return !(a > b);
}

template < typename RT > inline
bool operator<=(const Root_of_2<RT> &a, const RT &b) 
{
  return !(a > b);
}

template < typename RT > inline
bool operator<=(const CGAL_int(RT) &a, const Root_of_2<RT> &b) 
{
  return !(a > b);
}

template < typename RT > inline
bool operator<=(const Root_of_2<RT> &a, const CGAL_int(RT) &b) 
{
  return !(a > b);
}

template < typename RT > inline
bool operator==(const Root_of_2<RT> &a, const Root_of_2<RT> &b) 
{
  return CGAL_NTS compare(a, b) == 0;
}

template < typename RT > inline
bool operator==(const typename Root_of_traits< RT >::RootOf_1 &a,
	   const Root_of_2<RT> &b) 
{
  return CGAL_NTS compare(a, b) == 0;
}

template < typename RT > inline
bool operator==(const Root_of_2<RT> &a,
	   const typename Root_of_traits< RT >::RootOf_1 &b) 
{
  return CGAL_NTS compare(a, b) == 0;
}

template < typename RT > inline
bool operator==(const RT &a, const Root_of_2<RT> &b) 
{
  return CGAL_NTS compare(a, b) == 0;
}

template < typename RT > inline
bool operator==(const Root_of_2<RT> &a, const RT &b) 
{
  return CGAL_NTS compare(a, b) == 0;
}

template < typename RT > inline
bool operator==(const CGAL_int(RT) &a, const Root_of_2<RT> &b) 
{
  return CGAL_NTS compare(a, b) == 0;
}

template < typename RT > inline
bool operator==(const Root_of_2<RT> &a, const CGAL_int(RT) &b) 
{
  return CGAL_NTS compare(a, b) == 0;
}

template < typename RT > inline
bool operator!=(const Root_of_2<RT> &a, const Root_of_2<RT> &b) 
{
  return !(a == b);
}

template < typename RT > inline
bool operator!=(const typename Root_of_traits< RT >::RootOf_1 &a,
	   const Root_of_2<RT> &b) 
{
  return !(a == b);
}

template < typename RT > inline
bool operator!=(const Root_of_2<RT> &a,
	   const typename Root_of_traits< RT >::RootOf_1 &b) 
{
  return !(a == b);
}

template < typename RT > inline
bool operator!=(const RT &a, const Root_of_2<RT> &b) 
{
  return !(a == b);
}

template < typename RT > inline
bool operator!=(const Root_of_2<RT> &a, const RT &b) 
{
  return !(a == b);
}

template < typename RT > inline
bool operator!=(const CGAL_int(RT) &a, const Root_of_2<RT> &b) 
{
  return !(a == b);
}

template < typename RT > inline
bool operator!=(const Root_of_2<RT> &a, const CGAL_int(RT) &b) 
{
  return !(a == b);
}

// END OF COMPARISON OPERATORS

template < typename RT >
bool is_negative(const Root_of_2<RT> &a) 
{
  return sign(a) == NEGATIVE;
}

template < typename RT >
bool is_positive(const Root_of_2<RT> &a) 
{
  return sign(a) == POSITIVE;
}

template < typename RT >
bool is_zero(const Root_of_2<RT> &a) 
{
  return sign(a) == ZERO;
}


template < typename RT >
Root_of_2<RT>
square(const Root_of_2<RT> &a)
{

  if(a.is_rational()) {
    return Root_of_2<RT>(CGAL_NTS square(a[1]), CGAL_NTS square(a[2]));
  }

  if(is_zero(a[1])) {
    return Root_of_2<RT>(-a[0], a[2]);
  }

  CGAL_assertion(is_valid(a));

  // It's easy to get the explicit formulas for the square of the two roots.
  // Then it's easy to compute their sum and their product, which gives the
  // coefficients of the polynomial (X^2 - Sum X + Product).
 return Root_of_2<RT> ( CGAL_NTS square(a[2]),
                        2 * a[2] * a[0] - CGAL_NTS square(a[1]),
                        CGAL_NTS square(a[0]),
                        (a.is_smaller() ^ (a[1]>0)),
                        a.is_known_delta_not_zero());
}

template < typename RT >
Root_of_2<RT> inverse(const Root_of_2<RT> &a) 
{
  CGAL_assertion(is_valid(a));
  return a.inverse();
}

template < typename RT >
Root_of_2<RT> make_sqrt(const RT& r) 
{
  CGAL_assertion(r >= 0); 
  return Root_of_2<RT> (1,0,-r,false);
}

template < typename RT >
Root_of_2<RT> make_sqrt(const typename Root_of_traits< RT >::RootOf_1 r) 
{
  CGAL_assertion(r >= 0); 
  CGAL::Rational_traits< typename Root_of_traits< RT >::RootOf_1 > rat;
  return Root_of_2<RT> (rat.denominator(r),0,-rat.numerator(r),false);
}

// Mixed operators with Root_of_2 and RT/FT.
// Specializations of Binary_operator_result.
// Note : T1 can be different from T2 because of quotient types...
template < typename T1, typename T2 >
struct Binary_operator_result <Root_of_2<T1>, Root_of_2<T2> >;
// {
    // typedef void  type;
// };

template < typename T1, typename T2 >
struct Binary_operator_result <T1, Root_of_2<T2> > 
{
    typedef Root_of_2<T2>  type;
};

template < typename T1, typename T2 >
struct Binary_operator_result <Root_of_2<T1>, T2> 
{
    typedef Root_of_2<T1>  type;
};

template < typename RT >
struct Binary_operator_result <Root_of_2<RT>, typename Root_of_traits<RT>::RootOf_1 > 
{
    typedef Root_of_2<RT>  type;
};

template < typename RT >
struct Binary_operator_result <typename Root_of_traits<RT>::RootOf_1, Root_of_2<RT> > 
{
    typedef Root_of_2<RT>  type;
};








template < typename RT >
Root_of_2<RT>
operator-(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1& b)
{
  typedef typename Root_of_traits< RT >::RootOf_1  RootOf_1;
  typedef CGAL::Rational_traits< RootOf_1 >        Rational;
  //RT should be the same as Rational::RT

  Rational r;
  const RT &b1 = r.denominator(b);
  const RT &b0 = r.numerator(b);
  
  CGAL_assertion(is_valid(a) && is_valid(b));

  if(a.is_rational()) {
    if(is_zero(b0)) return Root_of_2<RT>(a[1],a[2]);
    return Root_of_2<RT>(a[1]*b1 - a[2]*b0,a[2]*b1);
  }

  RT sqb1 = CGAL_NTS square(b1);
  RT sqb0 = CGAL_NTS square(b0);
  RT b0b1 = b0 * b1;

  return Root_of_2<RT>(a[2] * sqb1,
                       2*a[2]*b0b1 + a[1]*sqb1,
                       a[2]* sqb0 + a[1]*b0b1 + a[0]*sqb1,
                       a.is_smaller(),
                       a.is_known_delta_not_zero());
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
  typedef typename Root_of_traits< RT >::RootOf_1  RootOf_1;

  CGAL_assertion(is_valid(a) && is_valid(b));

  if(a.is_rational()) {
    if(is_zero(b)) return Root_of_2<RT>(a[1],a[2]);
    return Root_of_2<RT>(a[1] - (b*a[2]),a[2]);
  }

  // It's easy to see it using the formula (X^2 - sum X + prod).

  //RT p = a[0] + b * a[1] + a[2] * CGAL_NTS square(b);
  //RT s = a[1] + b * a[2] * 2;

  RT tmp = a[2] * b;
  RT s = a[1] + tmp;
  RT p = a[0] + s * b;
  s += tmp;

  return Root_of_2<RT>(a[2], s, p, a.is_smaller(), a.is_known_delta_not_zero());
}

template < typename RT > inline
Root_of_2<RT> operator-(const Root_of_2<RT> &a, const CGAL_int(RT)& b) 
{
  return (a-RT(b));
}

template < typename RT > inline
Root_of_2<RT> operator-(const RT &a, const Root_of_2<RT> &b) 
{
  return (-(b-a));
}

template < typename RT > inline
Root_of_2<RT> operator-(const CGAL_int(RT)& a, const Root_of_2<RT> &b) 
{
  return (-(b-RT(a)));
}

template < typename RT > inline
Root_of_2<RT> operator+(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1& b) 
{
  return a - typename Root_of_traits< RT >::RootOf_1(-b);
}

template < typename RT > inline
Root_of_2<RT> operator+(const typename Root_of_traits< RT >::RootOf_1 &a,
	  const Root_of_2<RT> &b) 
{
  return b - typename Root_of_traits< RT >::RootOf_1(-a);
}

template < typename RT > inline
Root_of_2<RT> operator+(const Root_of_2<RT> &a, const RT& b) 
{
  return a - RT(-b);
}

template < typename RT > inline
Root_of_2<RT> operator+(const Root_of_2<RT> &a, const CGAL_int(RT)& b) 
{
  return a - RT(-b);
}

template < typename RT > inline
Root_of_2<RT> operator+(const RT &a, const Root_of_2<RT> &b) 
{
  return b - RT(-a);
}

template < typename RT > inline
Root_of_2<RT> operator+(const CGAL_int(RT) &a, const Root_of_2<RT> &b) 
{
  return b - RT(-a);
}

template < typename RT >
Root_of_2<RT>
operator*(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1& b)
{
  typedef typename Root_of_traits< RT >::RootOf_1  RootOf_1;
  typedef CGAL::Rational_traits< RootOf_1 >        Rational;
  //RT should be the same as Rational::RT

  Rational r;
  const RT &b1 = r.denominator(b);
  const RT &b0 = r.numerator(b);
  CGAL_assertion(is_valid(a) && is_valid(b));

  if(is_zero(b0)) return Root_of_2<RT>();

  if(a.is_rational()) {
    return Root_of_2<RT>(a[1]*b0,a[2]*b1);
  }

  return Root_of_2<RT>(a[2] * CGAL_NTS square(b1), a[1] * b1 * b0,
                       a[0] * CGAL_NTS square(b0),
                       CGAL_NTS sign(b) < 0 ? !a.is_smaller()
                                            :  a.is_smaller(),
                       a.is_known_delta_not_zero());
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
  typedef typename Root_of_traits< RT >::RootOf_1  RootOf_1;

  CGAL_assertion(is_valid(a));

  if(is_zero(b)) return Root_of_2<RT>();

  if(a.is_rational()) {
    if(is_zero(b)) return Root_of_2<RT>();
    return Root_of_2<RT>(a[1]*b,a[2]);
  }

  return Root_of_2<RT>(a[2], a[1] * b, a[0] * CGAL_NTS square(b),
                         b < 0 ? !a.is_smaller() : a.is_smaller());
}

template < typename RT > inline
Root_of_2<RT> operator*(const Root_of_2<RT> &a, const CGAL_int(RT)& b) 
{
  return a * RT(b);
}

template < typename RT > inline
Root_of_2<RT> operator*(const RT &a, const Root_of_2<RT> &b) 
{
  return b * a;
}

template < typename RT > inline
Root_of_2<RT> operator*(const CGAL_int(RT) &a, const Root_of_2<RT> &b) 
{
  return b * RT(a);
}

template < typename RT >
Root_of_2<RT>
operator/(const Root_of_2<RT> &a, const RT& b)
{
  typedef typename Root_of_traits< RT >::RootOf_1  RootOf_1;

  CGAL_assertion(b != 0);
  CGAL_assertion(is_valid(a));

  if(a.is_rational()) {
    return Root_of_2<RT>(a[1],a[2]*b);
  } 

  return Root_of_2<RT>(a[2] * CGAL_NTS square(b), a[1] * b, a[0],
                         b < 0 ? !a.is_smaller() : a.is_smaller(),
                       a.is_known_delta_not_zero());
}

template < typename RT > inline
Root_of_2<RT> operator/(const Root_of_2<RT> &a, const CGAL_int(RT)& b) 
{
  return a / RT(b);
}

template < typename RT > inline
Root_of_2<RT> operator/(const RT& a, const Root_of_2<RT> &b) 
{
  return b.inverse() * a;
}

template < typename RT > inline
Root_of_2<RT> operator/(const CGAL_int(RT)& a, const Root_of_2<RT> &b) 
{
  return b.inverse() * RT(a);
}

template < typename RT >
Root_of_2<RT>
operator/(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1& b)
{
  typedef typename Root_of_traits< RT >::RootOf_1  RootOf_1;
  typedef CGAL::Rational_traits< RootOf_1 >        Rational;
  //RT should be the same as Rational::RT

  CGAL_assertion(b != 0);
  CGAL_assertion(is_valid(a) && is_valid(b));

  Rational r;
  RT b0 = r.denominator(b);
  RT b1 = r.numerator(b);

  if(a.is_rational()) {
    return Root_of_2<RT>(a[1]*b0,a[2]*b1);
  } 

  if (b1<0)
    b0 = -b0, b1 = -b1;
  CGAL_assertion(b1>0);

  return Root_of_2<RT>(a[2] * CGAL_NTS square(b1), a[1] * b1 * b0,
                       a[0] * CGAL_NTS square(b0),
                       CGAL_NTS sign(b) < 0 ? !a.is_smaller()
                                            :  a.is_smaller(),
                       a.is_known_delta_not_zero());
}

template < typename RT > inline
Root_of_2<RT> operator/(const typename Root_of_traits< RT >::RootOf_1& a, 
          const Root_of_2<RT> &b) 
{
  return b.inverse() * a;
}

template < typename RT >
double
to_double(const Root_of_2<RT> &x)
{
  typedef typename Root_of_traits< RT >::RootOf_1  RootOf_1;

  CGAL_assertion(is_valid(x));

  if(x.is_rational()) {    
    double n = CGAL::to_double(x[1]);
    double d = CGAL::to_double(x[2]);
    return n/d;
  }

  if(is_zero(x[0])) {
    if(is_negative(x[1])) {
      if(x.is_smaller()) return 0.0;
      double a = CGAL::to_double(x[2]);
      double b = CGAL::to_double(x[1]);  
      return -b/a;
    } 
    if(x.is_smaller()) {
      double a = CGAL::to_double(x[2]);
      double b = CGAL::to_double(x[1]);  
      return -b/a;
    } return 0.0;
  }

  if(is_zero(x[1])) {
    double a = CGAL::to_double(x[2]);
    double c = CGAL::to_double(x[0]);
    if(x.is_smaller())
      return -CGAL::sqrt(-c/a);
    else return CGAL::sqrt(-c/a);
  }

  double a = CGAL::to_double(x[2]);
  double b = CGAL::to_double(x[1]);
  double d = std::sqrt(CGAL_NTS to_double(x.discriminant()));

  CGAL_assertion(a > 0);
  if (x.is_smaller())
    d = -d;

  return (d-b)/(a*2);
}

template < typename RT >
std::pair<double, double>
to_interval(const Root_of_2<RT> &x)
{
  CGAL_assertion(is_valid(x));

  if(x.is_rational()) {
    if(is_zero(x[1])) return std::make_pair(0.0,0.0);
    Interval_nt<> a = to_interval(x[1]);
    Interval_nt<> b = to_interval(x[2]);
    return (a/b).pair();
  }

  if(is_zero(x[0])) {
    if(is_negative(x[1])) {
      if(x.is_smaller()) return std::make_pair(0.0,0.0);
      Interval_nt<> a = to_interval(x[2]);
      Interval_nt<> b = to_interval(x[1]);  
      return (-b/a).pair();
    } 
    if(x.is_smaller()) {
      Interval_nt<> a = to_interval(x[2]);
      Interval_nt<> b = to_interval(x[1]);  
      return (-b/a).pair();
    } return std::make_pair(0.0,0.0);
  }

  if(is_zero(x[1])) {
    Interval_nt<> a = to_interval(x[2]);
    Interval_nt<> c = to_interval(x[0]);
    if(x.is_smaller())
      return (-CGAL::sqrt(-c/a)).pair();
    else return (CGAL::sqrt(-c/a)).pair();
  }

  Interval_nt<> a = to_interval(x[2]);
  Interval_nt<> b = to_interval(x[1]);
  Interval_nt<> c = to_interval(x[0]);
  Interval_nt<> d = CGAL::sqrt(CGAL_NTS square(b) - 4*a*c);

  if (x.is_smaller())
    d = -d;

  return ((d-b)/(a*2)).pair();
}

template < typename RT >
std::ostream &
operator<<(std::ostream &os, const Root_of_2<RT> &r)
{
  if(r.is_rational()) {
    return os << r.is_rational() << r[1] << " " << r[2];
  } else { 
    return os << r.is_rational() << r[2] << " "
	      << r[1] << " "
	      << r[0] << " "
	      << r.is_smaller() << " " << r.is_known_delta_not_zero(); 
  }
}

template < typename RT >
std::istream &
operator>>(std::istream &is, Root_of_2<RT> &r)
{
  RT a,b,c;
  bool s, izd;
  bool rat;
  is >> rat;
  if(rat) {
    is >> a >> b;
    if(is)
      r = Root_of_2<RT>(a,b);
    return is;
  }
  is >> a >> b >> c >> s >> izd;
  if(is)
    r = Root_of_2<RT>(a,b,c,s,izd);  
  return is;
}


template < typename RT >
void
print(std::ostream &os, const Root_of_2<RT> &r)
{
  if(r.is_rational()) {
    os << "Root_of_2( (" 
       << r[2] << ") * X - ("
       << r[1] << ") )";
  } else {
    os << "Root_of_2( (" 
       << r[2] << ") * X^2 + ("
       << r[1] << ") * X + ("
       << r[0] << ") , "
       << (r.is_smaller() ? "smaller" : "larger") << " )";
  }
}

template < typename RT >
bool
is_valid(const Root_of_2<RT> &r)
{
  return r.is_valid();
}


template < class RT >
struct Number_type_traits < Root_of_2<RT> >
{
  typedef Tag_false      Has_gcd;
  typedef Tag_false      Has_division;
  typedef Tag_false      Has_sqrt;
  typedef Tag_false      Has_exact_sqrt;
  typedef Tag_false      Has_exact_division;
  typedef Tag_false      Has_exact_ring_operations;
};

template < typename RT >
inline
io_Operator
io_tag(const Root_of_2<RT>&)
{
  return io_Operator();
}


} // namespace CGAL

#undef CGAL_int
#undef CGAL_double

#endif // CGAL_ROOT_OF_2_H
