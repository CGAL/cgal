// Copyright (c) 1999-2007
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra, Sylvain Pion, Michael Hemmer

// The template class Quotient<NT> is based on the LEDA class
// leda_rational written by Stefan Naeher and Christian Uhrig.
// It is basically a templated version with restricted functionality
// of the version of rational in LEDA release 3.3.
// The modification was done by Stefan.Schirra@mpi-sb.mpg.de

// The include is done before the protect macro on purpose, because
// of a cyclic dependency.

#ifndef CGAL_QUOTIENT_H
#define CGAL_QUOTIENT_H

#include <utility>
#include <istream>

#include <CGAL/Interval_nt.h>
#include <CGAL/Kernel/mpl.h>

#include <boost/operators.hpp>

namespace CGAL {

#define CGAL_int(T)    typename First_if_different<int,    T>::Type
#define CGAL_double(T) typename First_if_different<double, T>::Type

// Simplify the quotient numerator/denominator.
// Currently the default template doesn't do anything.
// This function is not documented as a number type requirement for now.
template < typename NT >
inline void
simplify_quotient(NT &, NT &) {}

// This one should be replaced by some functor or tag.
// Meanwhile, the class is specialized for Gmpz, mpz_class, leda_integer.
template < typename NT >
struct Split_double
{
  void operator()(double d, NT &num, NT &den) const
  {
    num = NT(d);
    den = 1;
  }
};


template <class NT_>
class Quotient
  : boost::ordered_field_operators1< Quotient<NT_>
  , boost::ordered_field_operators2< Quotient<NT_>, NT_
  , boost::ordered_field_operators2< Quotient<NT_>, CGAL_int(NT_)
  , boost::ordered_field_operators2< Quotient<NT_>, CGAL_double(NT_)
    > > > >
{
 public:
  typedef NT_        NT;

  Quotient()
    : num(0), den(1) {}

  Quotient(const NT& n)
    : num(n), den(1) {}

  Quotient(const CGAL_double(NT) & n)
  { Split_double<NT>()(n, num, den); }

  Quotient(const CGAL_int(NT) & n)
    : num(n), den(NT(1)) {}

  template <class T>
  explicit Quotient(const T& n) : num(n), den(1) {}

  template <class T>
  Quotient(const Quotient<T>& n) : num(n.numerator()), den(n.denominator()) {}

  Quotient& operator=(const NT & n)
  {
    num = n;
    den = 1;
    return *this;
  }

  Quotient& operator=(const CGAL_double(NT) & n)
  {
    Split_double<NT>()(n, num, den);
    return *this;
  }

  Quotient& operator=(const CGAL_int(NT) & n)
  {
    num = n;
    den = 1;
    return *this;
  }

  template <class T1, class T2>
  Quotient(const T1& n, const T2& d) : num(n), den(d)
  { CGAL_precondition( d != 0 ); }

  Quotient<NT>& operator+= (const Quotient<NT>& r);
  Quotient<NT>& operator-= (const Quotient<NT>& r);
  Quotient<NT>& operator*= (const Quotient<NT>& r);
  Quotient<NT>& operator/= (const Quotient<NT>& r);
  Quotient<NT>& operator+= (const NT& r);
  Quotient<NT>& operator-= (const NT& r);
  Quotient<NT>& operator*= (const NT& r);
  Quotient<NT>& operator/= (const NT& r);
  Quotient<NT>& operator+= (const CGAL_int(NT)& r);
  Quotient<NT>& operator-= (const CGAL_int(NT)& r);
  Quotient<NT>& operator*= (const CGAL_int(NT)& r);
  Quotient<NT>& operator/= (const CGAL_int(NT)& r);
  Quotient<NT>& operator+= (const CGAL_double(NT)& r);
  Quotient<NT>& operator-= (const CGAL_double(NT)& r);
  Quotient<NT>& operator*= (const CGAL_double(NT)& r);
  Quotient<NT>& operator/= (const CGAL_double(NT)& r);

  friend bool operator==(const Quotient& x, const Quotient& y)
  { return x.num * y.den == x.den * y.num; }
  friend bool operator==(const Quotient& x, const NT& y)
  { return x.den * y == x.num; }
  friend inline bool operator==(const Quotient& x, const CGAL_int(NT) & y)
  { return x.den * y == x.num; }
  friend inline bool operator==(const Quotient& x, const CGAL_double(NT) & y)
  { return x.den * y == x.num; } // Uh?

  Quotient<NT>&    normalize();

  const NT&   numerator()   const { return num; }
  const NT&   denominator() const { return den; }

  void swap(Quotient &q)
  {
    using std::swap;
    swap(num, q.num);
    swap(den, q.den);
  }

#ifdef CGAL_ROOT_OF_2_ENABLE_HISTOGRAM_OF_NUMBER_OF_DIGIT_ON_THE_COMPLEX_CONSTRUCTOR
  int tam() const { return (std::max)(num.tam(), den.tam()); }
#endif

 public:
  NT   num;
  NT   den;
};

template <class NT>
inline
void swap(Quotient<NT> &p, Quotient<NT> &q)
{
  p.swap(q);
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::normalize()
{
  if (num == den)
  {
      num = den = 1;
      return *this;
  }
  if (-num == den)
  {
      num = -1;
      den = 1;
      return *this;
  }
  NT ggt = CGAL_NTS gcd(num, den);
  if (ggt != 1 )
  {
      num = CGAL::integral_division(num, ggt);
      den = CGAL::integral_division(den, ggt);
  }
  return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator+= (const Quotient<NT>& r)
{
    num = num * r.den + r.num * den;
    den *= r.den;
    simplify_quotient(num, den);
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator-= (const Quotient<NT>& r)
{
    num = num * r.den - r.num * den;
    den *= r.den;
    simplify_quotient(num, den);
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator*= (const Quotient<NT>& r)
{
    num *= r.num;
    den *= r.den;
    simplify_quotient(num, den);
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator/= (const Quotient<NT>& r)
{
    CGAL_precondition( r.num != 0 );
    num *= r.den;
    den *= r.num;
    simplify_quotient(num, den);
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator+= (const NT& r)
{
    num += r * den;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator-= (const NT& r)
{
    num -= r * den;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator*= (const NT& r)
{
    num *= r;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator/= (const NT& r)
{
    CGAL_precondition( r != 0 );
    den *= r;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator+= (const CGAL_int(NT)& r)
{
    num += r * den;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator-= (const CGAL_int(NT)& r)
{
    num -= r * den;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator*= (const CGAL_int(NT)& r)
{
    num *= r;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator/= (const CGAL_int(NT)& r)
{
    CGAL_precondition( r != 0 );
    den *= r;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator+= (const CGAL_double(NT)& r)
{
  //num += r * den;
  NT r_num, r_den;
  Split_double<NT>()(r,r_num,r_den);
  num = num*r_den + r_num*den;
  den *=r_den;
  return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator-= (const CGAL_double(NT)& r)
{
  //num -= r * den;
  NT r_num, r_den;
  Split_double<NT>()(r,r_num,r_den);
  num =  num*r_den - r_num*den;
  den *= r_den;
  return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator*= (const CGAL_double(NT)& r)
{
  // num *= r;

  NT r_num, r_den;
  Split_double<NT>()(r,r_num,r_den);
  num *= r_num;
  den *= r_den;
  return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator/= (const CGAL_double(NT)& r)
{
  CGAL_precondition( r != 0 );
  NT r_num, r_den;
  Split_double<NT>()(r,r_num,r_den);
  num *= r_den;
  den *= r_num;
  return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Comparison_result
quotient_cmp(const Quotient<NT>& x, const Quotient<NT>& y)
{
    // No assumptions on the sign of  den  are made

    // code assumes that SMALLER == - 1;
    CGAL_precondition( SMALLER == static_cast<Comparison_result>(-1) );

    int xsign = CGAL_NTS sign(x.num) * CGAL_NTS sign(x.den) ;
    int ysign = CGAL_NTS sign(y.num) * CGAL_NTS sign(y.den) ;
    if (xsign == 0) return static_cast<Comparison_result>(-ysign);
    if (ysign == 0) return static_cast<Comparison_result>(xsign);
    // now (x != 0) && (y != 0)
    int diff = xsign - ysign;
    if (diff == 0)
    {
        int msign = CGAL_NTS sign(x.den) * CGAL_NTS sign(y.den);
        NT leftop  = NT(x.num * y.den * msign);
        NT rightop = NT(y.num * x.den * msign);
        return CGAL_NTS compare(leftop, rightop);
    }
    else
    {
        return (xsign < ysign) ? SMALLER : LARGER;
    }
}


template <class NT>
std::ostream&
operator<<(std::ostream& s, const Quotient<NT>& r)
{
   return s << r.numerator() << '/' << r.denominator();
}

template <class NT>
std::istream&
operator>>(std::istream& in, Quotient<NT>& r)
{
  /* format  num/den  or simply  num  */

  NT num,den=1;
  in >> num;
  if(!in) return in;
  std::istream::sentry s(in); // skip whitespace
  if(in.peek()!='/'){
          if(!in.good()){
                  in.clear(std::ios_base::eofbit);
                  // unlikely to be some other reason?
          }
  } else {
          char c;
          in.get(c); // remove the '/'
          in >> den;
          if(!in) return in;
  }
  r=Quotient<NT>(num,den);
  return in;
}

template< class NT >
inline
Quotient<NT>
operator+( const Quotient<NT>& x ) {
  return Quotient<NT>(x);
}

template <class NT>
inline
Quotient<NT>
operator-(const Quotient<NT>& x)
{ return Quotient<NT>(-x.num,x.den); }


template <class NT>
CGAL_MEDIUM_INLINE
NT
quotient_truncation(const Quotient<NT>& r)
{ return (r.num / r.den); }




template <class NT>
CGAL_MEDIUM_INLINE
bool
operator<(const Quotient<NT>& x, const Quotient<NT>& y)
{
  return quotient_cmp(x,y) == SMALLER;
}

template <class NT>
CGAL_MEDIUM_INLINE
bool
operator<(const Quotient<NT>& x, const NT& y)
{
  return quotient_cmp(x,Quotient<NT>(y)) == SMALLER;
}

template <class NT>
CGAL_MEDIUM_INLINE
bool
operator<(const Quotient<NT>& x, const CGAL_int(NT)& y)
{
  return quotient_cmp(x,Quotient<NT>(y)) == SMALLER;
}

template <class NT>
CGAL_MEDIUM_INLINE
bool
operator<(const Quotient<NT>& x, const CGAL_double(NT)& y)
{
  return quotient_cmp(x,Quotient<NT>(y)) == SMALLER;
}


template <class NT>
inline
bool
operator>(const Quotient<NT>& x, const NT& y)
{ return quotient_cmp(x,Quotient<NT>(y)) == LARGER; }

template <class NT>
inline
bool
operator>(const Quotient<NT>& x, const CGAL_int(NT)& y)
{ return quotient_cmp(x, Quotient<NT>(y)) == LARGER; }

template <class NT>
inline
bool
operator>(const Quotient<NT>& x, const CGAL_double(NT)& y)
{ return quotient_cmp(x, Quotient<NT>(y)) == LARGER; }


template< class NT >
class Is_valid< Quotient<NT> >
  : public CGAL::cpp98::unary_function< Quotient<NT>, bool > {
  public :
    bool operator()( const Quotient<NT>& x ) const {
      return is_valid(x.num) && is_valid(x.den);
    }
};


template <class NT>
inline
const NT&
denominator(const Quotient<NT>& q)
{ return q.den ; }

template <class NT>
inline
const NT&
numerator(const Quotient<NT>& q)
{ return q.num ; }

// The min/max are functions are needed since LEDA defines template
// min/max functions which clash with std::min/max with ADL.
template <class NT>
inline
const Quotient<NT>&
min
BOOST_PREVENT_MACRO_SUBSTITUTION
(const Quotient<NT>& p, const Quotient<NT>& q)
{
  return (std::min)(p, q);
}
template <class NT>
inline
const Quotient<NT>&
max
BOOST_PREVENT_MACRO_SUBSTITUTION
(const Quotient<NT>& p, const Quotient<NT>& q)
{
  return (std::max)(p, q);
}

/*
template <class NT>
NT
gcd(const NT&, const NT&)
{ return NT(1); }
*/

#undef CGAL_double
#undef CGAL_int

//
// Algebraic structure traits
//
namespace INTERN_QUOTIENT {
  template< class NT, class Sqrt_functor >
  class Sqrt_selector {
    public:
      class Sqrt
        : public CGAL::cpp98::unary_function< NT, NT > {
        public:
          NT operator()( const NT& x ) const {
            CGAL_precondition(x > 0);
            return NT(CGAL_NTS sqrt(x.numerator()*x.denominator()),
                      x.denominator());
          }
      };
  };

  template< class NT >
  class Sqrt_selector< NT, Null_functor > {
    public:
      typedef Null_functor Sqrt;
  };

// TODO: Algebraic_category could be Field_with_sqrt_tag, if NT
//       is INEXACT (because Sqrt can be inexact) and has a Sqrt-functor.
template<class NT> class Algebraic_structure_traits_quotient_base;

template< class NT > class Algebraic_structure_traits_quotient_base< Quotient<NT> >
  : public Algebraic_structure_traits_base< Quotient<NT>, Field_tag >  {

public:
    typedef Quotient<NT> Type;

    typedef typename Algebraic_structure_traits<NT>::Is_exact        Is_exact;
    typedef Tag_false Is_numerical_sensitive;



    class Is_square
        : public CGAL::cpp98::binary_function< Quotient<NT>, Quotient<NT>&, bool > {
    public:
        bool operator()( Quotient<NT> x, Quotient<NT>& y ) const {
            NT x_num, x_den, y_num, y_den;
            x.normalize();
            x_num = x.numerator();
            x_den = x.denominator();

            typename Algebraic_structure_traits<NT>::Is_square is_square;
            bool num_is_square = is_square(x_num,y_num);
            bool den_is_square = is_square(x_den,y_den);
            y= Quotient<NT>(y_num,y_den);
            return num_is_square && den_is_square;
        }
        bool operator()(Quotient<NT> x) const {
            x.normalize();
            typename Algebraic_structure_traits<NT>::Is_square is_square;
            return is_square(x.numerator())&&is_square(x.denominator());
        }

    };

    typedef typename boost::mpl::if_c<
        !boost::is_same< typename Algebraic_structure_traits<NT>::Sqrt,
                         Null_functor >::value,
         typename INTERN_QUOTIENT::Sqrt_selector< Type,
                                                  Is_exact >::Sqrt,
         Null_functor
                            >::type Sqrt;

    class Simplify
      : public CGAL::cpp98::unary_function< Type&, void > {
      public:
        void operator()( Type& x) const {
            x.normalize();
        }
    };
};


template<class NT> class Real_embeddable_traits_quotient_base;
// Real embeddable traits
template < class NT > class Real_embeddable_traits_quotient_base< Quotient<NT> >
  : public INTERN_RET::Real_embeddable_traits_base< Quotient<NT>,
                  typename Real_embeddable_traits< NT >::Is_real_embeddable > {
  public:
    typedef Quotient<NT> Type;

    class Compare
      : public CGAL::cpp98::binary_function< Type, Type,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Type& x,
                                            const Type& y ) const {
          return quotient_cmp(x, y);
        }
    };

    class To_double
      : public CGAL::cpp98::unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
        // Original global function was marked with an TODO!!
          if (x.num == 0 )
            return 0;

          double nd = CGAL_NTS to_double( x.num );

          if (x.den == 1 )
            return nd;

          double dd = CGAL_NTS to_double( x.den );

          if ( CGAL_NTS is_finite( x.den ) && CGAL_NTS is_finite( x.num ) )
            return nd/dd;

          if ( CGAL_NTS abs(x.num) > CGAL_NTS abs(x.den) )
          {
              NT  nt_div = x.num / x.den;
              double divd = CGAL_NTS to_double(nt_div);
              if ( divd >= std::ldexp(1.0,53) )
              { return divd; }
          }
          if ( CGAL_NTS abs(x.num) < CGAL_NTS abs(x.den) )
          { return 1.0 / CGAL_NTS to_double( NT(1) / x ); }

          return nd/dd;
        }
    };

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
          Interval_nt<> quot =
                          Interval_nt<>(CGAL_NTS to_interval(x.numerator())) /
                          Interval_nt<>(CGAL_NTS to_interval(x.denominator()));
          return std::make_pair(quot.inf(), quot.sup());
        }
    };

    class Is_finite
      : public CGAL::cpp98::unary_function< Type, bool > {
      public:
        bool operator()( const Type& x ) const {
          return CGAL_NTS is_finite(x.num) && CGAL_NTS is_finite(x.den);
        }
    };
};
} // namespace INTERN_QUOTIENT

template< class NT > class Algebraic_structure_traits< Quotient<NT> >
    : public INTERN_QUOTIENT::Algebraic_structure_traits_quotient_base<
Quotient<NT> >{};

template< class NT > class Real_embeddable_traits< Quotient<NT> >
    : public INTERN_QUOTIENT::Real_embeddable_traits_quotient_base<
Quotient<NT> >{};


// self coercion
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF_TEM( Quotient<NT>, class NT)

// from int to Quotient
template <class NT>
struct Coercion_traits<typename First_if_different<int, NT>::Type,Quotient<NT> >
{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef Quotient<NT> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const Quotient<NT>& x)   const { return x;}
        Type operator()(
                const typename First_if_different<int, NT>::Type& x) const {
            return Type(x);}
    };
};
template <class NT>
struct Coercion_traits<Quotient<NT>,typename First_if_different<int, NT>::Type>
    :public Coercion_traits<typename First_if_different<int, NT>::Type,
Quotient<NT> >{};

// from double to Quotient
template <class NT>
struct Coercion_traits<typename First_if_different<double, NT>::Type,
Quotient<NT> >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef Quotient<NT> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const Quotient<NT>& x)   const { return x;}
        Type operator()(
                const typename First_if_different<double, NT>::Type& x) const {
            return Type(x);}
    };
};
template <class NT>
struct Coercion_traits<Quotient<NT>,
typename First_if_different<double, NT>::Type>
    :public Coercion_traits<typename First_if_different<double, NT>::Type,
Quotient<NT> >
{};

// from NT to Quotient
CGAL_DEFINE_COERCION_TRAITS_FROM_TO_TEM ( NT, Quotient<NT>, class NT)

/*! \ingroup NiX_Fraction_traits_spec
 *  \brief Specialization of Fraction_traits for Quotient<NT>
 */
template <class NT>
class Fraction_traits< Quotient<NT> > {
public:
    typedef Quotient<NT> Type;
    typedef ::CGAL::Tag_true Is_fraction;
    typedef NT Numerator_type;
    typedef Numerator_type Denominator_type;

    //TODO: check whether Numerator_type has a GCD.
    //will use Scalar_factor from Scalar_factor_traits (not implemented yet)
    //for more details see EXACUS:NumeriX/include/NiX/Scalar_factor_traits.h
    typedef typename Algebraic_structure_traits< Numerator_type >::Gcd Common_factor;

    class Decompose {
    public:
        typedef Type first_argument_type;
        typedef Numerator_type& second_argument_type;
        typedef Numerator_type& third_argument_type;
        void operator () (
                const Type& rat,
                Numerator_type& num,
                Numerator_type& den) {
            num = rat.numerator();
            den = rat.denominator();
        }
    };

    class Compose {
    public:
        typedef Numerator_type first_argument_type;
        typedef Numerator_type second_argument_type;
        typedef Type result_type;
        Type operator ()(
                const Numerator_type& num ,
                const Numerator_type& den ) {
            Type result(num, den);
            return result;
        }
    };
};

} //namespace CGAL

namespace Eigen {
  template<class> struct NumTraits;
  template<class NT> struct NumTraits<CGAL::Quotient<NT> >
  {
    typedef CGAL::Quotient<NT> Real;
    typedef CGAL::Quotient<NT> NonInteger;
    typedef CGAL::Quotient<NT> Nested;
    typedef CGAL::Quotient<NT> Literal;

    static inline Real epsilon() { return NumTraits<NT>::epsilon(); }
    static inline Real dummy_precision() { return NumTraits<NT>::dummy_precision(); }

    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = NumTraits<NT>::RequireInitialization,
      ReadCost = 2*NumTraits<NT>::ReadCost,
      AddCost = 150,
      MulCost = 100
    };
  };
}

#endif  // CGAL_QUOTIENT_H
