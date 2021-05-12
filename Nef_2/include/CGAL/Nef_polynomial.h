// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_NEF_POLYNOMIAL_H
#define CGAL_NEF_POLYNOMIAL_H

#include <CGAL/license/Nef_2.h>


#include <CGAL/Nef_2/Polynomial.h>
#include <cstddef>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 3
#include <CGAL/Nef_2/debug.h>
#include <vector>

#include <CGAL/Kernel/mpl.h>
#include <CGAL/tss.h>

#include <boost/operators.hpp>

namespace CGAL {

#define CGAL_int(T)    typename First_if_different<int,    T>::Type
#define CGAL_double(T) typename First_if_different<double, T>::Type

template <class NT>
class Nef_polynomial
  : boost::ordered_field_operators1< Nef_polynomial<NT>
  , boost::ordered_field_operators2< Nef_polynomial<NT>, int
  > >
  , public Nef::Polynomial<NT>
{
  typedef typename CGAL::Nef::Polynomial<NT>  Base;
  typedef typename Base::size_type       size_type;

 protected:
  Nef_polynomial(size_type s) : Base(s) {}

 public:
  Nef_polynomial() : Base() {}
  Nef_polynomial(const NT& a0) : Base(a0) {}
  Nef_polynomial(const NT& a0, const NT& a1) : Base(a0,a1) {}
  Nef_polynomial(const NT& a0, const NT& a1, const NT& a2) : Base(a0,a1,a2) {}

  template <class Fwd_iterator>
  Nef_polynomial(std::pair<Fwd_iterator, Fwd_iterator> poly) : Base(poly) {}

  Nef_polynomial(CGAL_double(NT) n) : Base(n) {}
  Nef_polynomial(CGAL_double(NT) n1, CGAL_double(NT) n2) : Base(n1, n2) {}
  Nef_polynomial(CGAL_int(NT) n) : Base(NT(n)) {}
  Nef_polynomial(CGAL_int(NT) n1, CGAL_int(NT) n2) : Base(n1,n2) {}

  Nef_polynomial(const Base& p) : Base(p) {}

  Base & polynomial() { return static_cast<Base&>(*this); }
  const Base & polynomial() const  { return static_cast<const Base&>(*this); }

    static NT& infi_maximal_value() {
      CGAL_STATIC_THREAD_LOCAL_VARIABLE(NT, R_, 1);
      return R_;
    }

  friend bool operator==(const Nef_polynomial<NT> &a, const Nef_polynomial<NT> &b)
  {
    return a.polynomial() == b.polynomial();
  }

  friend bool operator==(const Nef_polynomial<NT> &a, const NT& b)
  {
    return a.polynomial() == b;
  }

  friend bool operator==(const Nef_polynomial<NT> &a, int b)
  {
    return a.polynomial() == b;
  }

  friend bool operator<(const Nef_polynomial<NT> &a, const Nef_polynomial<NT> &b)
  {
    return a.polynomial() < b.polynomial();
  }

  friend bool operator<(const Nef_polynomial<NT> &a, const NT& b)
  {
    return a.polynomial() < b;
  }

  friend bool operator<(const Nef_polynomial<NT> &a, int b)
  {
    return a.polynomial() < b;
  }

  friend bool operator>(const Nef_polynomial<NT> &a, int b)
  {
    return a.polynomial() > b;
  }
};

template <class NT>
inline
Nef_polynomial<NT> operator+(const Nef_polynomial<NT> &a)
{
  return a;
}

template <class NT>
inline
Nef_polynomial<NT> operator-(const Nef_polynomial<NT> &a)
{
  return - a.polynomial();
}

#undef CGAL_double
#undef CGAL_int


// TODO: integral_division to get it an UniqueFactorizationDomain
// TODO: div / mod  for EuclideanRing
template <class NT> class Algebraic_structure_traits< Nef_polynomial<NT> >
    : public Algebraic_structure_traits_base
             < Nef_polynomial<NT>, CGAL::Integral_domain_without_division_tag>
{
    typedef Algebraic_structure_traits<NT> AST_NT;
public:
    typedef Nef_polynomial<NT> Type;
    typedef typename AST_NT::Is_exact            Is_exact;
    typedef Tag_false                            Is_numerical_sensitive;
    class Integral_division
        : public CGAL::cpp98::binary_function< Type, Type,
                                Type > {
    public:
        Type operator()( const Type& x,
                const Type& y ) const {
          Type result = x / y;
          CGAL_postcondition_msg(result * y == x, "exact_division failed\n");
          return result;
        }
    };

    class Gcd
      : public CGAL::cpp98::binary_function< Type, Type, Type > {
    public:
        Type operator()( const Type& x, const Type& y ) const {
            // By definition gcd(0,0) == 0
          if( x == Type(0) && y == Type(0) )
            return Type(0);

          return CGAL::Nef::gcd( x, y );
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };
};

template <class NT> class Real_embeddable_traits< Nef_polynomial<NT> >
  : public INTERN_RET::Real_embeddable_traits_base< Nef_polynomial<NT> , CGAL::Tag_true > {
  public:
    typedef Nef_polynomial<NT> Type;
    class Abs
        : public CGAL::cpp98::unary_function< Type, Type> {
    public:
        Type inline operator()( const Type& x ) const {
            return (CGAL::Nef::sign( x ) == CGAL::NEGATIVE)? -x : x;
        }
    };

    class Sgn
      : public CGAL::cpp98::unary_function< Type, CGAL::Sign > {
      public:
        CGAL::Sign inline operator()( const Type& x ) const {
            return CGAL::Nef::sign( x );
        }
    };

    class Compare
      : public CGAL::cpp98::binary_function< Type, Type,
                                CGAL::Comparison_result > {
      public:
        CGAL::Comparison_result inline operator()(
                const Type& x,
                const Type& y ) const {
            return (CGAL::Comparison_result) CGAL::Nef::sign( x - y );
        }
    };

    class To_double
      : public CGAL::cpp98::unary_function< Type, double > {
      public:
        double inline operator()( const Type& p ) const {
            return CGAL::to_double(
                    p.eval_at(Nef_polynomial<NT>::infi_maximal_value()));
        }
    };

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& p ) const {
            return CGAL::to_interval(p.eval_at(Nef_polynomial<NT>::infi_maximal_value()));
        }
    };
};

template <typename NT>
inline Nef_polynomial<NT> min BOOST_PREVENT_MACRO_SUBSTITUTION
(const Nef_polynomial<NT>& x,const Nef_polynomial<NT>& y){
  return (x<=y)?x:y;
}

template <typename NT>
inline Nef_polynomial<NT> max BOOST_PREVENT_MACRO_SUBSTITUTION
(const Nef_polynomial<NT>& x,const Nef_polynomial<NT>& y){
  return (x>=y)?x:y;
}

template <typename NT>
class Fraction_traits<Nef_polynomial<NT> > {
public:
    typedef Nef_polynomial<NT> Type;
    typedef Fraction_traits<NT> Base_traits;
    typedef typename Base_traits::Is_fraction Is_fraction;
    typedef CGAL::Nef_polynomial<typename Base_traits::Numerator_type>
      Numerator_type;
    typedef typename Base_traits::Denominator_type Denominator_type;
    //TODO:    typedef Base_traits::Common_factor Common_factor;
    class Decompose {
    public:
        typedef Type first_argument_type;
        typedef Numerator_type second_argument_type;
        typedef Denominator_type third_argument_type;
        void operator () (const first_argument_type& rat,
                          second_argument_type& num,
                          third_argument_type& den) {
          typename Base_traits::Decompose decompose;
          third_argument_type num0;
          third_argument_type num1;
          third_argument_type den1;
          third_argument_type den0;
          decompose(rat[0], num0, den0);
          if(rat.degree() > 0) {
            decompose(rat[1], num1, den1);
            // TODO            den = den1/gcd(den0, den1)*den0;
            den = den1*den0;
            num = Numerator_type(num0*den1, num1*den0);
          } else {
            den = den0;
            num = Numerator_type(num0);
          }
        }
    };
    class Compose {
    public:
        typedef Numerator_type first_argument_type;
        typedef Denominator_type second_argument_type;
        typedef Type result_type;
        result_type operator () (const first_argument_type& num,
                                 const second_argument_type& den) {
          typename Base_traits::Compose compose;
          if(num.degree() == 0)
            return result_type(compose(num[0],den));
          else
            return result_type(compose(num[0],den),
                               compose(num[1],den));
        }
    };
};


} //namespace CGAL

#endif  // CGAL_NEF_POLYNOMIAL_H
