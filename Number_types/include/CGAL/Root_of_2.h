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
// Author(s)     : Sylvain Pion, Monique Teillaud, Athanasios Kakargias, Pedro Machado

#ifndef CGAL_ROOT_OF_2_H
#define CGAL_ROOT_OF_2_H

#include <iostream>
#include <CGAL/number_type_basic.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/NT_converter.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/enum.h>
#include <CGAL/tags.h>
#include <CGAL/Number_types/internal_functions_comparison_root_of_2.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/assertions.h>
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
// - add +, -, *, / with 2 root_of_2 (when it is possible - same gamma)
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
// - operator<<()
// - print() ("pretty" printing)
// - make_root_of_2()
// - add sqrt() (when it's degree 1), or a make_sqrt<RT>(const RT &r) ?
// - inverse()
// TODO :
// - use Boost.Operators.
// - add subtraction/addition with a degree 2 Root_of of the same field ?
// - add constructor from Polynomial ?
//   There should be a proper separate class Polynomial.
// - in compare_roots, we evaluate the polynomial at some FT, or at some
//   root of degree 1 polynomials.  It would be nice to have a separate
//   polynomial class which performed this task (and others)...
// - overloaded versions of make_root_of_2<>() for Lazy_exact_nt<> and others.

template <class RT> struct Root_of_traits;

template < typename RT_ >
class Root_of_2 {

  // the value is the root of P(X) = C2.X^2 + C1.X + C0,
  // and C2 > 0.

public:

  typedef RT_ RT;
  typedef typename Root_of_traits<RT>::RootOf_1 FT;

private:

  FT  _alpha, _beta, _gamma;
  bool _rational;

public:

#ifdef CGAL_ROOT_OF_2_ENABLE_HISTOGRAM_OF_NUMBER_OF_DIGIT_ON_THE_COMPLEX_CONSTRUCTOR
  static int max_num_digit;
  static int histogram[10000];
#endif

  Root_of_2()
    : _alpha(0), _rational(true)
  {
    CGAL_assertion(is_valid());
  }

  Root_of_2(const RT& c0)
    : _alpha(c0), _rational(true)
  {
    CGAL_assertion(is_valid());
  }

  Root_of_2(const typename First_if_different<int, RT>::Type & c0)
    : _alpha(RT(c0)), _rational(true)
  {
    CGAL_assertion(is_valid());
  }

  Root_of_2(const typename First_if_different<FT, RT>::Type & c0)
    : _alpha(c0), _rational(true)
  {
    CGAL_assertion(is_valid());
  }

  Root_of_2(const RT& a, const RT& b) {
    CGAL_assertion( b != 0 );
    _alpha = FT(a,b);
    _rational = true;
    CGAL_assertion(is_valid());
  }

  Root_of_2(const RT& a, const RT& b, const RT& c, const bool s)
  {
    if ( a != 0 ) {
      _alpha = FT(-b,2*a);
      _gamma = CGAL_NTS square(alpha()) - FT(c,a);
      if(CGAL_NTS is_zero(gamma())) {
	_rational = true;
      } else {
	_beta = (s ? -1 : 1);
	_rational = false;
      }
    }
    else {
      CGAL_assertion( b != 0 );
      _rational = true;
      _alpha = FT(-c,b);
      _beta = 0;
      _gamma = 0;
    }
    CGAL_assertion(is_valid());
  }

  Root_of_2(const typename First_if_different<FT, RT>::Type & c0,
            const typename First_if_different<FT, RT>::Type & c1,
            const typename First_if_different<FT, RT>::Type & c2)
  {
    if(CGAL_NTS is_zero(c1) || CGAL_NTS is_zero(c2)) {
      _alpha = c0;
      _rational = true;

#ifdef CGAL_ROOT_OF_2_ENABLE_HISTOGRAM_OF_NUMBER_OF_DIGIT_ON_THE_COMPLEX_CONSTRUCTOR
      int n_a = c0.tam();
      if(max_num_digit < n_a) max_num_digit = n_a;
      histogram[n_a]++;
#endif

    } else {
      _alpha = c0;
      _beta = c1;
      _gamma = c2;
      _rational = false;

#ifdef CGAL_ROOT_OF_2_ENABLE_HISTOGRAM_OF_NUMBER_OF_DIGIT_ON_THE_COMPLEX_CONSTRUCTOR
      int n_a = c0.tam();
      int n_b = c1.tam();
      int n_c = c2.tam();
      if(max_num_digit < n_a) max_num_digit = n_a;
      if(max_num_digit < n_b) max_num_digit = n_b;
      if(max_num_digit < n_c) max_num_digit = n_c;
      histogram[n_a]++;
      histogram[n_b]++;
      histogram[n_c]++;
#endif

    }
    CGAL_assertion(is_valid());
  }

  template <typename RT2>
  Root_of_2(const Root_of_2<RT2>& r)
    : _alpha(r.alpha()), _beta(r.beta()), _gamma(r.gamma()), _rational(r.is_rational())
  {
  }

  Root_of_2 operator-() const
  {
    if(is_rational()) return Root_of_2(-alpha());
    return Root_of_2 (-alpha(), -beta(), gamma());
  }

  bool is_valid() const
  {
    if(!is_rational()) {
      return gamma() >= 0;
    } return true;
  }

  bool is_rational() const
  {
    return _rational;
  }

  Root_of_2 inverse() const
  {
    CGAL_assertion(!(CGAL_NTS is_zero(alpha()) && (CGAL_NTS is_zero(beta()) || CGAL_NTS is_zero(gamma())))); // root = 0,
    FT r = (CGAL_NTS square(alpha())) -  (CGAL_NTS square(beta()))*gamma();
    CGAL_assertion(!(CGAL_NTS is_zero(r)
                   && (CGAL_NTS sign(beta()) != CGAL_NTS sign(alpha()))));
// ex. 6 - 2 sqrt(9)
    if(CGAL_NTS is_zero(r)) return Root_of_2(1 / (2 * alpha()));
    else return Root_of_2(alpha()/r, -beta()/r, gamma());
  }

  Root_of_2 conjugate() const
  {
    if(is_rational()) return Root_of_2(alpha());
    return Root_of_2(alpha(),-beta(),gamma());
  }

  const FT& alpha() const
  {
    return _alpha;
  }

  const FT& beta() const
  {
    return _beta;
  }

  const FT& gamma() const
  {
    return _gamma;
  }

  bool is_smaller() const
  {
    return beta() <= 0;
  }

  // The following functions deal with the internal polynomial.
  // Probably they should move to a separate polynomial class.

  RT operator[](int i) const
  {
    typedef Rational_traits< FT > Rational;
    CGAL_assertion((i>=0) & (i<3));
    Rational r;
    const RT r1 = r.numerator(alpha());
    const RT d1 = r.denominator(alpha());
    const RT r2 = r.numerator(beta());
    const RT d2 = r.denominator(beta());
    const RT r3 = r.numerator(gamma());
    const RT d3 = r.denominator(gamma());
    if(i == 0) {
      return (CGAL_NTS square(d2)) * d3;
    }
    if(i == 1) {
      return -2 * (CGAL_NTS square(d2)) * d3 * r1;
    }
    // i == 2
    return ((CGAL_NTS square(d2)) * d3 * (CGAL_NTS square(r1))) -
           ((CGAL_NTS square(d1)) * r3 * (CGAL_NTS square(r2)));

  }

  RT discriminant() const
  {
    if(is_rational()) return RT(0);
    return CGAL_NTS square(operator[](1)) -
           4*(operator[](0))*(operator[](2));
  }

  template < typename T >
  T eval_at(const T& x) const
  {
    if(is_rational()) return x * (operator[](0)) - (operator[](1));
    if(CGAL_NTS is_zero(x)) return (operator[](2));
    const bool zeroC0 = CGAL_NTS is_zero((operator[](2)));
    const bool zeroC1 = CGAL_NTS is_zero((operator[](1)));
    if(zeroC0 && zeroC1) return x * (operator[](0));
    if(zeroC0) return x * ((operator[](1)) + x * (operator[](0)));
    if(zeroC1) return (x * x * (operator[](0))) + (operator[](2));
    return (operator[](2)) + x * ((operator[](1)) + x * (operator[](0)));
  }

  template < typename T >
  Sign sign_at(const T &x) const
  {
    // Maybe there is slightly more efficient.
    return CGAL_NTS sign(eval_at(x));
  }

}; // Root_of_2

// COERCION_TRAITS BEGIN

CGAL_DEFINE_COERCION_TRAITS_FOR_SELF_TEM(Root_of_2<RT>,class RT)

template <class RT>
struct Coercion_traits< RT , Root_of_2<RT> >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef Root_of_2<RT> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const Root_of_2<RT>& x)   const { return x;}
        Type operator()(const RT& x) const {
            return Type(x);}
    };
};

template <class RT>
struct Coercion_traits< CGAL_int(RT) , Root_of_2<RT> >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef Root_of_2<RT> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const Root_of_2<RT>& x)   const { return x;}
        Type operator()(CGAL_int(RT) x) const {
            return Type(x);}
    };
};

template <class RT>
struct Coercion_traits< typename Root_of_traits<RT>::RootOf_1 , Root_of_2<RT> >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef Root_of_2<RT> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const Root_of_2<RT>& x)   const { return x;}
        Type operator()(const RT& x) const {
            return Type(x);}
    };
};


template <class RT, class NTX >
struct Coercion_traits< Root_of_2<RT> , NTX  >
    :public Coercion_traits<NTX , Root_of_2<RT> >{};

// COERCION_TRAITS END

#ifdef CGAL_ROOT_OF_2_ENABLE_HISTOGRAM_OF_NUMBER_OF_DIGIT_ON_THE_COMPLEX_CONSTRUCTOR

template < typename RT_ >
int Root_of_2<RT_>::max_num_digit = 0;

template < typename RT_ >
int Root_of_2<RT_>::histogram[10000];

#endif


template < class NT1,class NT2 >
struct NT_converter < Root_of_2<NT1> , Root_of_2<NT2> >
  : public std::unary_function< NT1, NT2 >
{
    Root_of_2<NT2>
    operator()(const Root_of_2<NT1> &a) const
    {
      if(a.is_rational()) {
        return Root_of_2<NT2>(NT_converter<NT1,NT2>()(a.alpha()));
      } else {
        return Root_of_2<NT2>(NT_converter<NT1,NT2>()(a.alpha()),
                              NT_converter<NT1,NT2>()(a.beta()),
                              NT_converter<NT1,NT2>()(a.gamma()));
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

template <class RT>
class Algebraic_structure_traits<Root_of_2<RT> >
    :public Algebraic_structure_traits_base<Root_of_2<RT> , Null_tag >{
public:

    typedef Root_of_2<RT> Type;
    typedef typename Algebraic_structure_traits<RT>::Is_exact Is_exact;
    struct Square
        : public std::unary_function< Root_of_2<RT> , Root_of_2<RT> >{
        Root_of_2<RT> operator()(const Root_of_2<RT>& a){

            CGAL_assertion(is_valid(a));

            if(a.is_rational()) {
                return Root_of_2<RT>(CGAL_NTS square(a.alpha()));
            }

            // It's easy to get the explicit formulas for the square of the two roots.
            // Then it's easy to compute their sum and their product, which gives the
            // coefficients of the polynomial (X^2 - Sum X + Product).
            return Root_of_2<RT> ( CGAL_NTS square(a.alpha()) +
                    (CGAL_NTS square(a.beta())) * a.gamma(),
                    2 * a.alpha() * a.beta(),
                    a.gamma());
        }
    };
};


template<class RT>
class Real_embeddable_traits<Root_of_2<RT> >
    :public INTERN_RET::Real_embeddable_traits_base<Root_of_2<RT> , CGAL::Tag_true >{
    typedef Real_embeddable_traits<RT> RET_RT;
    typedef typename Root_of_traits<RT>::RootOf_1 Root_of_1;
public:
    typedef Root_of_2<RT> Type;
    typedef Tag_true Is_real_embeddable;

    class Abs
        : public std::unary_function< Type, Type >{
    public:
        Type operator()(const Type& x) const {
            return (x>=0)?x:-x;
        }
    };

    class Sgn
        : public std::unary_function< Type, ::CGAL::Sign >{
    public:
        ::CGAL::Sign operator()(const Type& a) const {
            const ::CGAL::Sign sign_alpha = CGAL_NTS sign(a.alpha());
            if (a.is_rational()) return (sign_alpha);
            // If alpha and beta have the same sign, return this sign.
            const ::CGAL::Sign sign_beta = CGAL_NTS sign (a.beta());
            if (sign_alpha == sign_beta) return (sign_alpha);
            if (sign_alpha == ZERO) return (sign_beta);

            // Compare the squared values of m_alpha and of m_beta*sqrt(m_gamma):
            const Comparison_result res = CGAL_NTS compare (CGAL_NTS square(a.alpha()),
                    CGAL_NTS square(a.beta()) * a.gamma());
            if (res == LARGER) return sign_alpha;
            else if (res == SMALLER) return sign_beta;
            else return ZERO;
        }
    };

    class Compare
        : public std::binary_function< Type,
                                  Type,
                                  Comparison_result >{
    public:
        Comparison_result operator()(
                const Type& a,
                const Type& b) const{
            typedef typename Root_of_traits< RT >::RootOf_1 FT;
            typedef typename First_if_different<FT, RT>::Type WhatType;
            typedef typename boost::is_same< WhatType, RT > do_not_filter;

            CGAL_assertion(is_valid(a) & is_valid(b));

            if (a.is_rational()) return (CGAL_NTS compare(a.alpha(), b));
            if (b.is_rational()) return (CGAL_NTS compare(a, b.alpha()));

            if(!do_not_filter::value) {
                Interval_nt<> ia = CGAL_NTS to_interval(a);
                Interval_nt<> ib = CGAL_NTS to_interval(b);
                if(ia.inf() > ib.sup()) return LARGER;
                if(ia.sup() < ib.inf()) return SMALLER;
            }

            // Perform the exact comparison:
            // Note that the comparison of (a1 + b1*sqrt(c1)) and (a2 + b2*sqrt(c2))
            // is equivalent to comparing (a1 - a2) and (b2*sqrt(c2) -  b1*sqrt(c1)).
            // We first determine the signs of these terms.
            const FT diff_alpha = a.alpha() - b.alpha();
            const ::CGAL::Sign sign_left = CGAL_NTS sign(diff_alpha);
            const FT a_sqr = a.beta()*a.beta()*a.gamma();
            const FT b_sqr = b.beta()*b.beta()*b.gamma();
            Comparison_result right_res = CGAL_NTS compare (b_sqr, a_sqr);
            ::CGAL::Sign sign_right = ZERO;

            if (right_res == LARGER)
                {
                    // Take the sign of b2:
                    sign_right = CGAL_NTS sign(b.beta());
                }
            else if (right_res == SMALLER)
                {
                    // Take the opposite sign of b1:
                    switch (CGAL_NTS sign (a.beta()))
                        {
                        case POSITIVE :
                            sign_right = NEGATIVE;
                            break;
                        case NEGATIVE:
                            sign_right = POSITIVE;
                            break;
                        case ZERO:
                            sign_right = ZERO;
                            break;
                        default:
                            // We should never reach here.
                            CGAL_error();
                        }
                }
            else
                {
                    // We take the sign of (b2*sqrt(c2) -  b1*sqrt(c1)), where both terms
                    // have the same absolute value. The sign is equal to the sign of b2,
                    // unless both terms have the same sign, so the whole expression is 0.
                    sign_right = CGAL_NTS sign (b.beta());
                    if (sign_right == CGAL_NTS sign (a.beta()))
                        sign_right = ZERO;
                }

            // Check whether on of the terms is zero. In this case, the comparsion
            // result is simpler:
            if (sign_left == ZERO)
                {
                    if (sign_right == POSITIVE)
                        return (SMALLER);
                    else if (sign_right == NEGATIVE)
                        return (LARGER);
                    else
                        return (EQUAL);
                }
            else if (sign_right == ZERO)
                {
                    if (sign_left == POSITIVE)
                        return (LARGER);
                    else if (sign_left == NEGATIVE)
                        return (SMALLER);
                    else
                        return (EQUAL);
                }

            // If the signs are not equal, we can determine the comparison result:
            if (sign_left != sign_right)
                {
                    if (sign_left == POSITIVE)
                        return (LARGER);
                    else
                        return (SMALLER);
                }

            // We now square both terms and look at the sign of the one-root number:
            //   ((a1 - a2)^2 - (b1^2*c1 + b2^2*c2)) + 2*b1*b2*sqrt(c1*c2)
            //
            // If both signs are negative, we should swap the comparsion result
            // we eventually compute.
            const FT A = diff_alpha*diff_alpha - (a_sqr + b_sqr);
            const FT B = 2 * a.beta() * b.beta();
            const FT C = a.gamma() * b.gamma();
            const ::CGAL::Sign sgn = CGAL_NTS sign(Root_of_2<RT>(A, B, C));
            const bool swap_res = (sign_left == NEGATIVE);

            if (sgn == POSITIVE)
                return (swap_res ? SMALLER : LARGER);
            else if (sgn == NEGATIVE)
                return (swap_res ? LARGER : SMALLER);
            else
                return (EQUAL);
        }

        Comparison_result
        inline
        operator()(
                const Type& a,
                const Root_of_1& b
        ) const{
            typedef typename Root_of_traits< RT >::RootOf_1 FT;
            typedef typename First_if_different<FT, RT>::Type WhatType;
            typedef typename boost::is_same< WhatType, RT > do_not_filter;

            CGAL_assertion(is_valid(a) & is_valid(b));

            if (a.is_rational()) return (CGAL_NTS compare (a.alpha(), b));

            if(!do_not_filter::value) {
                Interval_nt<> ia = CGAL_NTS to_interval(a);
                Interval_nt<> ib = CGAL_NTS to_interval(b);
                if(ia.inf() > ib.sup()) return LARGER;
                if(ia.sup() < ib.inf()) return SMALLER;
            }

            // Perform the exact comparison.
            const ::CGAL::Sign   sgn = CGAL_NTS  sign(a - b);
            if (sgn == POSITIVE) return (LARGER);
            else if (sgn == NEGATIVE) return (SMALLER);
            else return (EQUAL);
        }

        Comparison_result
        inline
        operator()(
                const Root_of_1& a,
                const Type& b
        ) const{ return opposite(this->operator()(b,a) ); }

        Comparison_result
        inline
        operator()(
                const Type& a,
                const RT& b
        ) const{
            typedef typename Root_of_traits< RT >::RootOf_1 FT;
            typedef typename First_if_different<FT, RT>::Type WhatType;
            typedef typename boost::is_same< WhatType, RT > do_not_filter;

            CGAL_assertion(is_valid(a) & is_valid(b));

            if (a.is_rational()) return (CGAL_NTS compare (a.alpha(), b));

            if(!do_not_filter::value) {
                Interval_nt<> ia = CGAL_NTS to_interval(a);
                Interval_nt<> ib = CGAL_NTS to_interval(b);
                if(ia.inf() > ib.sup()) return LARGER;
                if(ia.sup() < ib.inf()) return SMALLER;
            }

            // Perform the exact comparison.
            const ::CGAL::Sign   sgn = CGAL_NTS sign(a - b);
            if (sgn == POSITIVE) return (LARGER);
            else if (sgn == NEGATIVE) return (SMALLER);
            else return (EQUAL);
        }

        inline
        Comparison_result
        operator()(const RT &a, const Root_of_2<RT> &b)
        {
            return opposite(this->operator()(b, a));
        }

        inline
        Comparison_result
        operator()( CGAL_int(RT)  a, const Root_of_2<RT> &b)
        {
            return this->operator()(RT(a),b);
        }

        inline
        Comparison_result
        operator()(const Root_of_2<RT> &a, CGAL_int(RT) b)
        {
            return this->operator()(a,RT(b));
        }
    };

    class To_double
        : public std::unary_function< Type, double >{
    public:
        double operator()(const Type& x) const {
            if (x.is_rational()) {
			        return (CGAL_NTS to_double(x.alpha()));
			      }
            return CGAL_NTS to_double(x.alpha()) +
                   CGAL_NTS to_double(x.beta()) *
                   (std::sqrt)(CGAL_NTS to_double(x.gamma()));
        }
    };

    class To_interval
        : public std::unary_function< Type, std::pair< double, double > >{
    public:
        std::pair<double,double> operator()(const Type& x) const {

            if(x.is_rational()) return CGAL_NTS to_interval(x.alpha());

            const Interval_nt<true>   alpha_in
                = CGAL_NTS to_interval(x.alpha());
            const Interval_nt<true>   beta_in
                = CGAL_NTS to_interval(x.beta());
            const Interval_nt<true>   gamma_in
                = CGAL_NTS to_interval(x.gamma());
            const Interval_nt<true>&  x_in = alpha_in +
                (beta_in * CGAL_NTS sqrt(gamma_in));
            return x_in.pair();
        }
    };
};

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
Root_of_2<RT> inverse(const Root_of_2<RT> &a)
{
  CGAL_assertion(is_valid(a));
  return a.inverse();
}

template < typename RT >
Root_of_2<RT> make_sqrt(const RT& r)
{
  CGAL_assertion(r >= 0);
  if(CGAL_NTS is_zero(r)) return Root_of_2<RT>();
  return Root_of_2<RT>(0,1,r);
}

template < typename RT >
Root_of_2<RT> make_sqrt(const typename Root_of_traits< RT >::RootOf_1& r)
{
  CGAL_assertion(r >= 0);
  if(CGAL_NTS is_zero(r)) return Root_of_2<RT>();
  return Root_of_2<RT>(0,1,r);
}


template < typename RT >
Root_of_2<RT>
operator-(const Root_of_2<RT> &a,
	  const typename Root_of_traits< RT >::RootOf_1& b)
{
  typedef typename Root_of_traits< RT >::RootOf_1  RootOf_1;
  typedef Rational_traits< RootOf_1 >        Rational;
  //RT should be the same as Rational::RT

  CGAL_assertion(is_valid(a) & is_valid(b));

  if(a.is_rational()) {
    return Root_of_2<RT>(a.alpha() - b);
  }

  return Root_of_2<RT>(a.alpha() - b, a.beta(), a.gamma());
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
  typedef Rational_traits< RootOf_1 >        Rational;
  //RT should be the same as Rational::RT

  CGAL_assertion(is_valid(a) & is_valid(b));

  if(a.is_rational()) {
    return Root_of_2<RT>(a.alpha() - b);
  }

  return Root_of_2<RT>(a.alpha() - b, a.beta(), a.gamma());
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
  typedef Rational_traits< RootOf_1 >        Rational;
  //RT should be the same as Rational::RT

  CGAL_assertion(is_valid(a) & is_valid(b));

  if(CGAL_NTS is_zero(b)) return Root_of_2<RT>();

  if(a.is_rational()) {
    return Root_of_2<RT>(a.alpha() * b);
  }

  return Root_of_2<RT>(a.alpha() * b,
                       a.beta() * b,
                       a.gamma());
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
  typedef Rational_traits< RootOf_1 >        Rational;
  //RT should be the same as Rational::RT

  CGAL_assertion(is_valid(a) & is_valid(b));

  if(CGAL_NTS is_zero(b)) return Root_of_2<RT>();

  if(a.is_rational()) {
    return Root_of_2<RT>(a.alpha() * b);
  }

  return Root_of_2<RT>(a.alpha() * b,
                       a.beta() * b,
                       a.gamma());
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
    return Root_of_2<RT>(a.alpha() / b);
  }

  return Root_of_2<RT>(a.alpha()/b,
                       a.beta()/b,
                       a.gamma());
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

  CGAL_assertion(b != 0);
  CGAL_assertion(is_valid(a));

  if(a.is_rational()) {
    return Root_of_2<RT>(a.alpha() / b);
  }

  return Root_of_2<RT>(a.alpha()/b,
                       a.beta()/b,
                       a.gamma());
}

template < typename RT > inline
Root_of_2<RT> operator/(const typename Root_of_traits< RT >::RootOf_1& a,
          const Root_of_2<RT> &b)
{
  return b.inverse() * a;
}

template < typename RT >
Root_of_2<RT>
operator-(const Root_of_2<RT> &a,
	  const Root_of_2<RT> &b)
{

  CGAL_assertion(is_valid(a));
  CGAL_assertion(is_valid(b));
  CGAL_assertion((a.is_rational() || b.is_rational()) || (a.gamma() == b.gamma()));

  if(a.is_rational() && b.is_rational()) {
    return Root_of_2<RT>(a.alpha() - b.alpha());
  }
  if(a.is_rational()) return a.alpha() - b;
  if(b.is_rational()) return a - b.alpha();

  if(a.beta() == b.beta()) {
    return Root_of_2<RT>(a.alpha() - b.alpha());
  }

  return Root_of_2<RT>(a.alpha() - b.alpha(),
                       a.beta() - b.beta(),
                       a.gamma());
}

template < typename RT > inline
Root_of_2<RT> operator+(const Root_of_2<RT> &a, const Root_of_2<RT> &b)
{
  return b - (-a);
}

template < typename RT >
Root_of_2<RT>
operator*(const Root_of_2<RT> &a,
	  const Root_of_2<RT> &b)
{

  CGAL_assertion(is_valid(a));
  CGAL_assertion(is_valid(b));
  CGAL_assertion((a.is_rational() || b.is_rational()) || (a.gamma() == b.gamma()));

  if(a.is_rational() && b.is_rational()) {
    return Root_of_2<RT>(a.alpha() * b.alpha());
  }

  if(a.is_rational()) {
    if(CGAL_NTS is_zero(a.alpha())) return Root_of_2<RT>();
    return Root_of_2<RT>(b.alpha() * a.alpha(), b.beta() * a.alpha(), b.gamma());
  }

  if(b.is_rational()) {
    if(CGAL_NTS is_zero(b.alpha())) return Root_of_2<RT>();
    return Root_of_2<RT>(a.alpha() * b.alpha(), a.beta() * b.alpha(), a.gamma());
  }

  return Root_of_2<RT>(b.beta() * a.beta() * a.gamma() + a.alpha() * b.alpha(),
                       a.alpha() * b.beta() + a.beta() * b.alpha(),
                       a.gamma());
}

template < typename RT > inline
Root_of_2<RT> operator/(const Root_of_2<RT> &a, const Root_of_2<RT> &b)
{
  return b.inverse() * a;
}

template < typename RT >
double
to_double(const Root_of_2<RT> &x)
{
  if (x.is_rational()) {
    return (CGAL_NTS to_double(x.alpha()));
  }
  return CGAL_NTS to_double(x.alpha()) +
         CGAL_NTS to_double(x.beta()) *
         (std::sqrt)(CGAL_NTS to_double(x.gamma()));
}

template < typename RT >
std::ostream &
operator<<(std::ostream &os, const Root_of_2<RT> &r)
{
  if(r.is_rational()) {
    return os << r.is_rational() << " " << r.alpha();
  } else {
    return os << r.is_rational() << " " << r.alpha() << " "
	      << r.beta() << " "
	      << r.gamma();
  }
}

template < typename RT >
std::istream &
operator>>(std::istream &is, Root_of_2<RT> &r)
{
  typedef typename Root_of_traits< RT >::RootOf_1  FT;
  FT a,b,c;
  bool rat;
  is >> rat;
  if(rat) {
    is >> a;
    if(is) r = Root_of_2<RT>(a);
    return is;
  }
  is >> a >> b >> c;
  if(is) r = Root_of_2<RT>(a,b,c);
  return is;
}


template < typename RT >
void
print(std::ostream &os, const Root_of_2<RT> &r)
{
  if(r.is_rational()) {
    os << "(" << r.alpha() << ")";
  } else {
    os << "(" << r.alpha() << " + " << r.beta() <<
          "*sqrt(" << r.gamma() << ")"<< ")";
  }
}

template < typename RT >
class Is_valid<Root_of_2<RT> >: public std::unary_function<Root_of_2<RT> , bool>{
public:
    bool operator()(const Root_of_2<RT> &r)
    {
        return r.is_valid();
    }
};

template <class NT>
inline const Root_of_2<NT>& min BOOST_PREVENT_MACRO_SUBSTITUTION
(const Root_of_2<NT>& p, const Root_of_2<NT>& q){
  return (std::min)(p, q);
}
template <class NT> 
inline const Root_of_2<NT>& max BOOST_PREVENT_MACRO_SUBSTITUTION
(const Root_of_2<NT>& p, const Root_of_2<NT>& q){
  return (std::max)(p, q);
}

} // namespace CGAL

#undef CGAL_int
#undef CGAL_double

#endif // CGAL_ROOT_OF_2_H
