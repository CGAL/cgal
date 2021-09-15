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

#if defined(CGAL_USE_CPP_INT) || !defined(CGAL_DO_NOT_RUN_TESTME)
#include <CGAL/boost_mp.h>
// TODO: The two below will be removed when we move simplify_quotient() to boost_mp,
// where it is supposed to be.
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#endif

#include <boost/operators.hpp>
#include <boost/type_index.hpp>

namespace CGAL {

#define CGAL_int(T)    typename First_if_different<int,    T>::Type
#define CGAL_double(T) typename First_if_different<double, T>::Type

// Simplify the quotient numerator/denominator.
// Currently the default template doesn't do anything.
// This function is not documented as a number type requirement for now.
template < typename NT >
inline void
simplify_quotient(NT & a, NT & b) { }

#if defined(CGAL_USE_CPP_INT) || !defined(CGAL_DO_NOT_RUN_TESTME)

// WARNINGS:
// https://cgal.geometryfactory.com/~cgaltest/test_suite/TESTRESULTS/CGAL-5.4-Ic-5937/Convex_decomposition_3/TestReport_cgaltest_Ubuntu-gcc7.gz
// https://cgal.geometryfactory.com/~cgaltest/test_suite/TESTRESULTS/CGAL-5.4-Ic-5937/Convex_decomposition_3/TestReport_cgaltest_arm64-apple_bigsur_clang1200-release-64bits.gz
// See boost includes above. Should make it work, I guess.
inline void
simplify_quotient(boost::multiprecision::cpp_int & a, boost::multiprecision::cpp_int & b) {

  // TODO:
  // - move it to the boost_mp.h
  // - can we use gcd only sometimes to save time?
  const boost::multiprecision::cpp_int r = boost::multiprecision::gcd(a, b);
  a /= r;
  b /= r;
}

#endif // CPP_INT

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

#if defined(CGAL_USE_CPP_INT) || !defined(CGAL_DO_NOT_RUN_TESTME)

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator+= (const Quotient<NT>& r)
{
#if true

    num = num * r.den + r.num * den;
    den *= r.den;
    simplify_quotient(num, den);
    return *this;

#else

    // taken from https://github.com/boostorg/multiprecision/blob/rational_adaptor_standalone/include/boost/multiprecision/rational_adaptor.hpp#L486
    // This has the BSL, but this code won't survive
    typedef NT Backend;

    NT result_num, result_denom;
    Backend  gcd, t1, t2, t3, t4;
   //
   // Begin by getting the gcd of the 2 denominators:
   //
   // ERROR: no matching function for call to gcd! We skip this implementation for the moment!
    gcd = boost::multiprecision::gcd(denominator(), b.denominator());
   //
   // Do we have gcd > 1:
   //
   // WARNING: logical not is only applied to the left hand side of comparison and
   // ERROR: no match for operator!
   // We skip this implementation for the moment!
   if (! gcd == 1)
   {
      //
      // Scale the denominators by gcd, and put the results in t1 and t2:
      //
     t1 = CGAL::integral_division( b.denominator(), gcd);
     t2 = CGAL::integral_division( denominator(), gcd);
      //
      // multiply the numerators by the scale denominators and put the results in t3, t4:
      //
      t3 =  numerator() * t1;
      t4 =  b.numerator() * t2;
      //
      // Add them up:
      //
      t3 += t4;
      //
      // Get the gcd of gcd and our numerator (t3):
      //
      // ERROR: no matching function for call to gcd! We skip this implementation for the moment!
      t4 = boost::multiprecision::gcd(t3, gcd);
      if (t4 == 1)
      {
         result_num = t3;
         result_denom = t1 * denominator();
      }
      else
      {
         //
         // Uncommon case where gcd is not 1, divide the numerator
         // and the denominator terms by the new gcd.  Note we perform division
         // on the existing gcd value as this is the smallest of the 3 denominator
         // terms we'll be multiplying together, so there's a good chance it's a
         // single limb value already:
         //
  result_num= CGAL::integral_division(t3, t4);
        t3 = CGAL::integral_division( gcd, t4);
         t4 =  t1 * t2;
         result_denom =  t4 * t3;
      }
   }
   else
   {
      //
      // Most common case (approx 60%) where gcd is one:
      //
      t1 = numerator() * b.denominator();
      t2 =  denominator() *b.numerator();

      result_num  =  t1 + t2;

      result_denom = denominator() * b.denominator();

   }

  *this = Quotient<NT>(result_num,result_denom);
  return *this;
#endif

}

#else // default impl

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

#endif // default impl

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

        bool are_correct_bounds( const double l, const double u, const Type& x ) const {

          if (
            CGAL::abs(l) == std::numeric_limits<double>::infinity() ||
            CGAL::abs(u) == std::numeric_limits<double>::infinity() ||
            CGAL::abs(l) == 0.0 ||
            CGAL::abs(u) == 0.0) {
            return true;
          }
          CGAL_assertion(CGAL::abs(l) != std::numeric_limits<double>::infinity());
          CGAL_assertion(CGAL::abs(u) != std::numeric_limits<double>::infinity());
          CGAL_assertion(CGAL::abs(l) != 0.0);
          CGAL_assertion(CGAL::abs(u) != 0.0);

          const Type lb(l), ub(u);
          CGAL_assertion(lb <= x);
          CGAL_assertion(ub >= x);
          const bool are_bounds_respected = (lb <= x && x <= ub);

          const double inf = std::numeric_limits<double>::infinity();
          CGAL_assertion(u == l || u == std::nextafter(l, +inf));
          bool are_bounds_tight = (u == l || u == std::nextafter(l, +inf));
          return are_bounds_respected && are_bounds_tight;
        }

        #if defined(CGAL_USE_CPP_INT) || !defined(CGAL_DO_NOT_RUN_TESTME)

        std::pair< Interval_nt<>, int64_t > get_interval_exp( NT& x ) const {

          CGAL_assertion(CGAL::is_positive(x));
          int64_t d = 0;
          double l = 0.0, u = 0.0;
          const int64_t n = static_cast<int64_t>(boost::multiprecision::msb(x)) + 1;
          const int64_t num_dbl_digits = std::numeric_limits<double>::digits;

          if (n > num_dbl_digits) {
            d = n - num_dbl_digits;
            x >>= d;
            const uint64_t xx = static_cast<uint64_t>(x);
            const uint64_t yy = xx + 1;
            CGAL_assertion(xx > 0 && yy > xx);
            l = static_cast<double>(xx);
            u = static_cast<double>(yy);
          } else {
            const uint64_t xx = static_cast<uint64_t>(x);
            CGAL_assertion(xx > 0);
            l = static_cast<double>(xx);
            u = l;
          }
          return std::make_pair( Interval_nt<>(l, u), d );
        }

        std::pair<double, double> get_interval_as_gmpzf( Type x ) const {

          CGAL_assertion_code(const Type input = x);
          double l = 0.0, u = 0.0;
          if (CGAL::is_zero(x.num)) { // return [0.0, 0.0]
            CGAL_assertion(are_correct_bounds(l, u, input));
            return std::make_pair(l, u);
          }
          CGAL_assertion(!CGAL::is_zero(x.num));
          CGAL_assertion(!CGAL::is_zero(x.den));

          // Handle signs.
          bool change_sign = false;
          const bool is_num_pos = CGAL::is_positive(x.num);
          const bool is_den_pos = CGAL::is_positive(x.den);
          if (!is_num_pos && !is_den_pos) {
            x.num = -x.num;
            x.den = -x.den;
          } else if (!is_num_pos && is_den_pos) {
            change_sign = true;
            x.num = -x.num;
          } else if (is_num_pos && !is_den_pos) {
            change_sign = true;
            x.den = -x.den;
          }
          CGAL_assertion(CGAL::is_positive(x.num) && CGAL::is_positive(x.den));

          const auto num = get_interval_exp(x.num);
          const auto den = get_interval_exp(x.den);

          const Interval_nt<> div = num.first / den.first;
          const int64_t e = num.second - den.second;
          std::tie(l, u) = ldexp(div, e).pair();

          if (change_sign) {
            const double t = l;
            l = -u;
            u = -t;
          }

          CGAL_assertion(are_correct_bounds(l, u, input));
          return std::make_pair(l, u);
        }

        std::pair<double, double> get_1ulp_interval( const int64_t shift, const NT& p ) const {

          std::cout << "- 1ulp interval: " << std::endl;

          CGAL_assertion(p >= 0);
          const uint64_t pp = static_cast<uint64_t>(p);
          const uint64_t qq = static_cast<uint64_t>(pp + 1);

          CGAL_assertion(pp >= 0);
          CGAL_assertion(qq > pp);
          const double pp_dbl = static_cast<double>(pp);
          const double qq_dbl = static_cast<double>(qq);

          std::cout << "pp_dbl: " << pp_dbl << std::endl;
          std::cout << "qq_dbl: " << qq_dbl << std::endl;

          const Interval_nt<> intv(pp_dbl, qq_dbl);
          return ldexp(intv, -shift).pair();
        }

        std::pair<double, double> get_0ulp_interval( const int64_t shift, const NT& p ) const {

          std::cout << "- 0ulp interval: " << std::endl;

          CGAL_assertion(p >= 0);
          const uint64_t pp = static_cast<uint64_t>(p);

          CGAL_assertion(pp >= 0);
          const double pp_dbl = static_cast<double>(pp);

          std::cout << "pp_dbl: " << pp_dbl << std::endl;
          std::cout << "qq_dbl: " << pp_dbl << std::endl;

          const Interval_nt<> intv(pp_dbl, pp_dbl);
          return ldexp(intv, -shift).pair();
        }

        std::pair<double, double> get_interval_as_boost( Type x ) const {

          CGAL_assertion_code(const Type input = x);
          double l = 0.0, u = 0.0;
          if (CGAL::is_zero(x.num)) { // return [0.0, 0.0]
            CGAL_assertion(are_correct_bounds(l, u, input));
            return std::make_pair(l, u);
          }
          CGAL_assertion(!CGAL::is_zero(x.num));
          CGAL_assertion(!CGAL::is_zero(x.den));

          // Handle signs.
          bool change_sign = false;
          const bool is_num_pos = CGAL::is_positive(x.num);
          const bool is_den_pos = CGAL::is_positive(x.den);
          if (!is_num_pos && !is_den_pos) {
            x.num = -x.num;
            x.den = -x.den;
          } else if (!is_num_pos && is_den_pos) {
            change_sign = true;
            x.num = -x.num;
          } else if (is_num_pos && !is_den_pos) {
            change_sign = true;
            x.den = -x.den;
          }
          CGAL_assertion(CGAL::is_positive(x.num) && CGAL::is_positive(x.den));

          const int64_t num_dbl_digits = std::numeric_limits<double>::digits - 1;
          const int64_t msb_num = static_cast<int64_t>(boost::multiprecision::msb(x.num));
          const int64_t msb_den = static_cast<int64_t>(boost::multiprecision::msb(x.den));
          const int64_t msb_diff = msb_num - msb_den;
          int64_t shift = num_dbl_digits - msb_diff;

          std::cout << "msb  num: " << msb_num  << std::endl;
          std::cout << "msb  den: " << msb_den  << std::endl;
          std::cout << "msb diff: " << msb_diff << std::endl;
          std::cout << "   shift: " << shift    << std::endl;

          std::cout << "x num: " << x.num << std::endl;
          std::cout << "x den: " << x.den << std::endl;

          if (shift > 0) {
            CGAL_assertion(msb_diff < num_dbl_digits);
            x.num <<= shift;
          } else if (shift < 0) {
            CGAL_assertion(msb_diff > num_dbl_digits);
            x.den <<= boost::multiprecision::detail::unsigned_abs(shift);
          }
          CGAL_assertion(num_dbl_digits ==
            static_cast<int64_t>(boost::multiprecision::msb(x.num)) -
            static_cast<int64_t>(boost::multiprecision::msb(x.den)));

          NT p, r;
          boost::multiprecision::divide_qr(x.num, x.den, p, r);
          CGAL_assertion_code(const int64_t p_bits =
            static_cast<int64_t>(boost::multiprecision::msb(p)));

          std::cout << "x num shifted: " << x.num << std::endl;
          std::cout << "x den shifted: " << x.den << std::endl;

          std::cout << "p bits: " << boost::multiprecision::msb(p) << std::endl;

          std::cout << "p: " << p << std::endl;
          std::cout << "r: " << r << std::endl;

          if (r == 0) {
            std::cout << "- case r = 0" << std::endl;
            std::tie(l, u) = get_0ulp_interval(shift, p);
          } else {
            if (p_bits == num_dbl_digits - 1) {

              std::cout << "- case r > 0 && p_bits = 51" << std::endl;
              r <<= 1;
              const int cmp = r.compare(x.den);
              if (cmp > 0) {

                p <<= 1;
                ++shift;
                ++p;
                std::tie(l, u) = get_1ulp_interval(shift, p);

              } else if ( ((cmp == 0) && (p & 1u))) {
                CGAL_assertion_msg(false, "TODO: THIS CASE1!");
              } else {
                CGAL_assertion_msg(false, "TODO: THIS CASE2!");
              }

            } else {
              std::cout << "- case r > 0 && p_bits = 52" << std::endl;
              std::tie(l, u) = get_1ulp_interval(shift, p);
            }
          }

          if (change_sign) {
            const double t = l;
            l = -u;
            u = -t;
          }

          std::cout << "          l: " << l << std::endl;
          std::cout << "          u: " << u << std::endl;
          const double inf = std::numeric_limits<double>::infinity();
          std::cout << "nextafter l: " << nextafter(l, +inf) << std::endl;
          std::cout << std::endl;

          CGAL_assertion(are_correct_bounds(l, u, input));
          return std::make_pair(l, u);
        }

        std::pair<double, double> interval_from_cpp_rational( const Type& x ) const {

          // Seems fast enough because this conversion happens
          // only a few times during the run, at least for NEF.
          boost::multiprecision::cpp_rational rat;
          CGAL_assertion(!CGAL::is_zero(x.den));
          if (CGAL::is_negative(x.den)) {
            rat = boost::multiprecision::cpp_rational(-x.num, -x.den);
          } else {
            CGAL_assertion(CGAL::is_positive(x.den));
            rat = boost::multiprecision::cpp_rational( x.num,  x.den);
          }

          double l, u;
          std::tie(l, u) = to_interval(rat);
          const double inf = std::numeric_limits<double>::infinity();

          if (l == +inf) {
            l = std::numeric_limits<double>::max();
            CGAL_assertion(u == +inf);
          } else if (u == -inf) {
            u = std::numeric_limits<double>::lowest();
            CGAL_assertion(l == -inf);
          }

          CGAL_assertion(are_correct_bounds(l, u, x));
          return std::make_pair(l, u);
        }

        std::pair<double, double> get_interval_using_cpp_rational( const Type& x ) const {

          const double inf = std::numeric_limits<double>::infinity();

          const Interval_nt<> xn = Interval_nt<>(CGAL_NTS to_interval(x.num));
          if (CGAL::abs(xn.inf()) == inf || CGAL::abs(xn.sup()) == inf) {
            return interval_from_cpp_rational(x);
          }
          CGAL_assertion(CGAL::abs(xn.inf()) != inf && CGAL::abs(xn.sup()) != inf);

          const Interval_nt<> xd = Interval_nt<>(CGAL_NTS to_interval(x.den));
          if (CGAL::abs(xd.inf()) == inf || CGAL::abs(xd.sup()) == inf) {
            return interval_from_cpp_rational(x);
          }
          CGAL_assertion(CGAL::abs(xd.inf()) != inf && CGAL::abs(xd.sup()) != inf);

          const Interval_nt<> quot = xn / xd;
          CGAL_assertion(are_correct_bounds(quot.inf(), quot.sup(), x));
          return std::make_pair(quot.inf(), quot.sup());
        }

        #endif // CPP_INT

        std::pair<double, double> operator()( const Type& x ) const {

        #if defined(CGAL_USE_CPP_INT) || !defined(CGAL_DO_NOT_RUN_TESTME)

          // Option 1.
          // Seems to be less precise and we rarely end up with an interval [d,d]
          // even for numbers, which are exactly representable as double.
          // Otherwise, it is quite similar to the results of the option 3.
          // return get_interval_as_gmpzf(x);

          // Option 2.
          // Works slightly better than the first one.
          // It always returns a correct interval, but it is sometimes less tight
          // than the one from the option 3.
          return get_interval_as_boost(x);

          // Option 3.
          // This seems to give the tightest intervals.
          // return get_interval_using_cpp_rational(x);

        #else // master version

          #if defined(CGAL_USE_TO_INTERVAL_WITH_BOOST)
          return std::make_pair(0.0, 0.0);
          #else
          const Interval_nt<> quot =
            Interval_nt<>(CGAL_NTS to_interval(x.numerator())) /
            Interval_nt<>(CGAL_NTS to_interval(x.denominator()));
          return std::make_pair(quot.inf(), quot.sup());
          #endif

        #endif // master version
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
