// Copyright (c) 1999-2005  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// Author(s)     : Stefan Schirra, Sylvain Pion

// The template class Quotient<NT> is based on the LEDA class
// leda_rational written by Stefan Naeher and Christian Uhrig.
// It is basically a templated version with restricted functionality
// of the version of rational in LEDA release 3.3.
// The modification was done by Stefan.Schirra@mpi-sb.mpg.de

// The include is done before the protect macro on purpose, because
// of a cyclic dependency.
#include <CGAL/basic.h>

#ifndef CGAL_QUOTIENT_H
#define CGAL_QUOTIENT_H

#include <utility>

#ifndef CGAL_CFG_NO_LOCALE
#  include <locale>
#else
#  include <cctype>
#endif

#include <CGAL/Quotient_fwd.h>

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/utils.h>

#include <CGAL/Interval_nt.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Kernel/mpl.h>

#include <boost/operators.hpp>
#include <boost/mpl/if.hpp>

#include <CGAL/Root_of_traits.h>
#include <CGAL/make_root_of_2.h>

#include <CGAL/number_utils.h>

#include <CGAL/functional_base.h> // Unary_function, Binary_function

CGAL_BEGIN_NAMESPACE

#define CGAL_int(T)    typename First_if_different<int,    T>::Type
#define CGAL_double(T) typename First_if_different<double, T>::Type

// Simplify the quotient numerator/denominator.
// Currently the default template doesn't do anything.
// This function is not documented as a number type requirement for now.
template < typename NT >
inline void
simplify_quotient(NT &, NT &) {}

template <class NT_>
class Quotient
  : boost::ordered_field_operators1< Quotient<NT_>
  , boost::ordered_field_operators2< Quotient<NT_>, NT_
  , boost::ordered_field_operators2< Quotient<NT_>, CGAL_int(NT_)
    > > >
{
 public:
  typedef NT_        NT;
  typedef Tag_false  Has_gcd;
  typedef Tag_true   Has_division;
  typedef typename Number_type_traits<NT_>::Has_sqrt  Has_sqrt;

  typedef Tag_true   Has_exact_division;
  typedef typename Number_type_traits<NT_>::Has_exact_sqrt Has_exact_sqrt;
  typedef typename Number_type_traits<NT_>::Has_exact_ring_operations
  Has_exact_ring_operations;

  Quotient()
    : num(0), den(1) {}

  Quotient(const NT& n)
    : num(n), den(1) {}

  Quotient(const CGAL_double(NT) & n)
    : num(n), den(1) {}

  Quotient(const CGAL_int(NT) & n)
    : num(n), den(1) {}

  template <class T>
  explicit Quotient(const T& n) : num(n), den(1) {}

  template <class T>
  Quotient(const Quotient<T>& n) : num(n.numerator()), den(n.denominator()) {}

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
  int tam() const { return std::max(num.tam(), den.tam()); }
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
  NT ggt = gcd(num, den);
  if (ggt != 1 )
  {
      num /= ggt;
      den /= ggt;
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
        NT leftop  = x.num * y.den * msign;
        NT rightop = y.num * x.den * msign;
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
   return s << r.numerator() << "/" << r.denominator();
}

template <class NT>
std::istream&
operator>>(std::istream& in, Quotient<NT>& r)
{
  /* format  num/den  or simply  num  */

  char c = 0;

#ifndef CGAL_CFG_NO_LOCALE
  while (in.get(c) && std::isspace(c, std::locale::classic() ));
#else
  while (in.get(c) && CGAL_CLIB_STD::isspace(c));
#endif // CGAL_CFG_NO_LOCALE
  if ( !in ) return in;
  in.putback(c);

  NT num;
  NT den(1);
  in >> num;

#ifndef CGAL_CFG_NO_LOCALE
  while (in.get(c) && std::isspace(c, std::locale::classic() ));
#else
  while (in.get(c) && CGAL_CLIB_STD::isspace(c));
#endif // CGAL_CFG_NO_LOCALE
  if (( in ) && ( c == '/'))
  {
#ifndef CGAL_CFG_NO_LOCALE
      while (in.get(c) && std::isspace(c, std::locale::classic() ));
#else
      while (in.get(c) && CGAL_CLIB_STD::isspace(c));
#endif // CGAL_CFG_NO_LOCALE
      CGAL_assertion( in );
      in.putback(c);
      in >> den;
  }
  else
  {
      in.putback(c);
      if ( in.eof() ) in.clear();
  }
  r = Quotient<NT>( num, den);
  return in;
}

template <class NT>
inline
io_Operator
io_tag(const Quotient<NT>&)
{ return io_Operator(); }


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
operator==(const Quotient<NT>& x, const Quotient<NT>& y)
{ return x.num * y.den == x.den * y.num; }

template <class NT>
CGAL_MEDIUM_INLINE
bool
operator==(const Quotient<NT>& x, const NT& y)
{ return x.den * y == x.num; }

template <class NT>
inline
bool
operator==(const Quotient<NT>& x, const CGAL_int(NT) & y)
{ return x.den * y == x.num; }



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
inline
bool
operator>(const Quotient<NT>& x, const NT& y)
{ return quotient_cmp(x,Quotient<NT>(y)) == LARGER; }

template <class NT>
inline
bool
operator>(const Quotient<NT>& x, const CGAL_int(NT)& y)
{ return quotient_cmp(x, Quotient<NT>(y)) == LARGER; }


template< class NT >
class Is_valid< Quotient<NT> > 
  : public Unary_function< Quotient<NT>, bool > {
  public :
    bool operator()( const Quotient<NT>& x ) {
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

template < class NT >
Quotient<NT> exact_division(const Quotient<NT>& n, const Quotient<NT>& d) {
    return n/d;
}

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
        : public Unary_function< NT, NT > {
        public:
          NT operator()( const NT& x ) const {
            CGAL_precondition(x > 0);
            return NT(CGAL_NTS sqrt(x.numerator()*x.denominator()),
                      x.denominator());          
          }
      };
  };
  
  template< class NT >
  class Sqrt_selector< NT, CGAL::Null_functor > {
    public:
      typedef CGAL::Null_functor Sqrt;
  };

// TODO: Algebraic_structure_tag could be Field_with_sqrt_tag, if NT
//       is INEXACT (because Sqrt can be inexact) and has a Sqrt-functor.
template<class NT> class Algebraic_structure_traits_quotient_base;

template< class NT > class Algebraic_structure_traits_quotient_base< Quotient<NT> >
  : public Algebraic_structure_traits_base< Quotient<NT>, CGAL::Field_tag >  {
public:
    typedef Quotient<NT> Algebraic_structure;  
 
    typedef typename Algebraic_structure_traits<NT>::Is_exact        Is_exact;
    

    
    class Is_square
        : public Binary_function< Quotient<NT>, Quotient<NT>&, bool > {
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
                         CGAL::Null_functor >::value,
         typename INTERN_QUOTIENT::Sqrt_selector< Algebraic_structure, 
                                                  Is_exact >::Sqrt,
         CGAL::Null_functor 
                            >::type Sqrt;

    class Simplify 
      : public Unary_function< Algebraic_structure&, void > {
      public:
        void operator()( Algebraic_structure& x) const {
            x.normalize();
        }
    };
};


template<class NT> class Real_embeddable_traits_quotient_base;
// Real embeddable traits
template < class NT > class Real_embeddable_traits_quotient_base< Quotient<NT> > 
  : public INTERN_RET::Real_embeddable_traits_base_selector< Quotient<NT>,
                  typename Real_embeddable_traits< NT >::Is_real_embeddable > {
  public:
    typedef Quotient<NT> Real_embeddable;  
               
    class Compare 
      : public Binary_function< Real_embeddable, Real_embeddable,
                                CGAL::Comparison_result > {
      public:
        CGAL::Comparison_result operator()( const Real_embeddable& x, 
                                            const Real_embeddable& y ) const {
          return quotient_cmp(x, y);
        }
    };
    
    class To_double 
      : public Unary_function< Real_embeddable, double > {
      public:
        double operator()( const Real_embeddable& x ) const {
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
              if ( divd >= CGAL_CLIB_STD::ldexp(1.0,53) )
              { return divd; }
          }
          if ( CGAL_NTS abs(x.num) < CGAL_NTS abs(x.den) )
          { return 1.0 / CGAL_NTS to_double( NT(1) / x ); }
        
          return nd/dd;
        }
    };
    
    class To_interval 
      : public Unary_function< Real_embeddable, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Real_embeddable& x ) const {
          Interval_nt<> quot = 
                          Interval_nt<>(CGAL_NTS to_interval(x.numerator())) /
                          Interval_nt<>(CGAL_NTS to_interval(x.denominator()));
          return std::make_pair(quot.inf(), quot.sup());
        }
    };
    
    class Is_finite
      : public Unary_function< Real_embeddable, bool > {
      public:
        bool operator()( const Real_embeddable& x ) {
          return CGAL_NTS is_finite(x.num) && CGAL_NTS is_finite(x.den);
        }
    };
};
} // namespace INTERN_QUOTIENT 

template< class NT > class Algebraic_structure_traits< Quotient<NT> >
    : public INTERN_QUOTIENT::Algebraic_structure_traits_quotient_base<Quotient<NT> >{};

template< class NT > class Real_embeddable_traits< Quotient<NT> >
    : public INTERN_QUOTIENT::Real_embeddable_traits_quotient_base<Quotient<NT> >{};


// fwd 
template <class A, class B> class Coercion_traits;
// self coercion
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF_TEM( Quotient<NT>, class NT);
// from int to Quotient

template <class NT>                                                      
struct Coercion_traits<typename First_if_different<int, NT>::Type,Quotient<NT> >{                                    
    typedef CGAL::Tag_true  Are_explicit_interoperable;             
    typedef CGAL::Tag_true  Are_implicit_interoperable;             
    typedef Quotient<NT> Coercion_type;                                       
    struct Cast{                                                    
        typedef Coercion_type result_type;                          
        Coercion_type operator()(const Quotient<NT>& x)   const { return x;}  
        Coercion_type operator()(const typename First_if_different<int, NT>::Type& x) const {             
            return Coercion_type(x);}                               
    };                                                              
};                                                                  
template <class NT>                                                      
struct Coercion_traits<Quotient<NT>,typename First_if_different<int, NT>::Type>
    :public Coercion_traits<typename First_if_different<int, NT>::Type,Quotient<NT> >
{};

// from double to Quotient
template <class NT>                                                      
struct Coercion_traits<typename First_if_different<double, NT>::Type,Quotient<NT> >{                                    
    typedef CGAL::Tag_true  Are_explicit_interoperable;             
    typedef CGAL::Tag_true  Are_implicit_interoperable;             
    typedef Quotient<NT> Coercion_type;                                       
    struct Cast{                                                    
        typedef Coercion_type result_type;                          
        Coercion_type operator()(const Quotient<NT>& x)   const { return x;}  
        Coercion_type operator()(const typename First_if_different<double, NT>::Type& x) const {             
            return Coercion_type(x);}                               
    };                                                              
};                                                                  
template <class NT>                                                      
struct Coercion_traits<Quotient<NT>,typename First_if_different<double, NT>::Type>
    :public Coercion_traits<typename First_if_different<double, NT>::Type,Quotient<NT> >
{};





// from NT to Quotient
CGAL_DEFINE_COERCION_TRAITS_FROM_TO_TEM ( NT, Quotient<NT>, class NT);




// Rational traits
template < class NT >
struct Rational_traits< Quotient<NT> >
{
  typedef NT RT;

  const RT & numerator   (const Quotient<NT>& r) const { return r.numerator(); }
  const RT & denominator (const Quotient<NT>& r) const { return r.denominator(); }
  
  Quotient<NT> make_rational(const RT & n, const RT & d) const
  { return Quotient<NT>(n, d); } 
  Quotient<NT> make_rational(const Quotient<NT> & n,
                             const Quotient<NT> & d) const
  { return n / d; } 
};

template < class NT >
inline
typename Root_of_traits< NT >::RootOf_2
make_root_of_2(const Quotient< NT > &a, const Quotient< NT > &b,
               const Quotient< NT > &c, bool d)
{
  return CGALi::make_root_of_2_rational< NT, Quotient< NT > >(a,b,c,d);
}

template < class NT >
inline
typename Root_of_traits< NT >::RootOf_2
make_root_of_2(const Quotient< NT > &a, const Quotient< NT > &b, const Quotient< NT > &c)
{
  return typename Root_of_traits< NT >::RootOf_2(a,b,c);
}

// CGAL::Quotient<NT> should be the same as Root_of_traits<NT>::RootOf_1
// i.e the default implementation.
template < class NT >
struct Root_of_traits< Quotient< NT > >
  : public Root_of_traits< NT > {};

CGAL_END_NAMESPACE

#endif  // CGAL_QUOTIENT_H
