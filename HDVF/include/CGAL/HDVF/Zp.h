// Copyright (c) 2025 LIS Marseille (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alexandra Bac <alexandra.bac@univ-amu.fr>


#ifndef CGAL_HDVF_ZP_H
#define CGAL_HDVF_ZP_H

#include <CGAL/license/HDVF.h>

#include <CGAL/number_type_basic.h>
#include <CGAL/number_utils.h>

#include <iostream>

namespace CGAL {
namespace Homological_discrete_vector_field {

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Zp` implements the concept `IntegralDomainWithoutDivision` with the ring \f$\mathbb Z/p\mathbb Z\f$ (which is a field when `p`is prime). This is a "lightweight" implementation which aims at providing fast operations and constructors.

 According to the value of `p`, users can chose the size of the representation used to store values (default size: `int`).

 \warning For \f$\mathbb Z/2\mathbb Z\f$, prefer the class `Z2` which is optimized.

 \cgalModels{IntegralDomainWithoutDivision}

 \tparam p a positive integer.
 \tparam T a type used for the inner storage of the values (default: `int`).
 \tparam IsPrime a boolean encoding weather `p` is a prime number or not.
 */

template <size_t p, typename T = unsigned int, bool IsPrime = true>
class Zp {
    T _i ;
public:

    /** \brief Constructor from a value */
    Zp(T i=0) : _i( (i>=0)?(i % T(p)):((i % T(p)) + T(p)) ) { }

    // Copy constructor
    Zp(const Zp& a) : _i(a._i) {}

    /** \brief Returns the value of p. */
    static T operator() () { return T(p); }
    
    /** \brief Tests if the element is 0. */
    bool is_zero() const { return _i == 0 ; }

    /** \brief Operator=. */
    Zp& operator= (const Zp& a) {
        _i = a._i;
        return *this;
    }
    
    /** \brief Unary operator+ */
    friend Zp operator+ (const Zp& a)
    {
        return Zp(a) ;
    }

    /** \brief Unary operator-. */
    friend Zp operator- (const Zp& a)
    {
        return Zp(T(p) - a._i) ;
    }

    /** \brief Operator+. */
    friend Zp     operator+ (const Zp& a, const Zp& b)
    {
        return Zp<p, T, IsPrime>(a._i + b._i) ;
    }

    /** \brief Operator-. */
    friend Zp     operator- (const Zp& a, const Zp& b) {
        if (a._i >= b._i)
            return Zp<p, T, IsPrime>(a._i - b._i) ;
        else
            return Zp<p, T, IsPrime>((T(p) - b._i) + a._i);
    }

    /** \brief Operator*. */
    friend Zp     operator* (const Zp& a, const Zp& b)
    {
        return Zp<p, T, IsPrime>(a._i * b._i) ;
    }

    /** \brief Operator/. */
    friend Zp     operator/ (const Zp& a, const Zp& b)
    {
        return Zp<p, T, IsPrime>(a._i / b._i) ;
    }

    /** \brief Operator+=. */
    Zp &     operator+= (const Zp& a)
    {
        _i += a._i;
        _i %= T(p) ;
        return *this ;
    }

    /** \brief Operator-=. */
    Zp &     operator-= (const Zp& a)
    {
        if (_i >= a._i)
            _i -= a._i ;
        else
        {
            _i += (T(p) - a._i) ;
        }
        return *this ;
    }

    /** \brief Operator*=. */
    Zp &     operator*= (const Zp& a)
    {
        _i *= a._i ;
        _i %= T(p) ;
        return *this ;
    }

    /** \brief Operator/=. */
    Zp &     operator/= (const Zp& a)
    {
        _i /= a._i ;
        _i %= T(p) ;
        return *this ;
    }

    /** \brief Operator==. */
    friend bool     operator== (const Zp& a, const Zp& b)
    {
        return (a._i == b._i) ;
    }

    /** \brief Operator!=. */
    friend bool     operator!= (const Zp& a, const Zp& b)
    {
        return (a._i != b._i);
    }

    /** \brief Absolute value. */
    friend Zp  abs(const Zp& a)
    {
        return Zp<p,T, IsPrime>(a) ;
    }

    /** \brief Operator<<. */
    friend std::ostream& operator<<(std::ostream& out, const Zp& a)
    {
        return (out << int(a._i)) ;
    }

    /** \brief Operator>>. */
    friend std::istream& operator>>(std::istream& in, Zp& a)
    {
        int tmp ;
        in >> tmp ;
        a = Zp(tmp) ;
        return (in) ;
    }
    
    /** \brief Returns the invertibility of the element. */
    bool is_invertible() {
        if (IsPrime)
            return (_i != 0);
        else
        {
            // TODO: optimize with static data
            return (std::gcd(_i,p) == 1);
        }
    }
};

} /* end namespace Homological_discrete_vector_field */

// Specialization for p being not a prime number
template <int p, typename T> class Algebraic_structure_traits< Homological_discrete_vector_field::Zp<p, T, false> >
  : public Algebraic_structure_traits_base< Homological_discrete_vector_field::Zp<p, T, false>, Integral_domain_without_division_tag >  {
  public:
      typedef Tag_true            Is_exact;
      typedef Tag_false           Is_numerical_sensitive;
      typedef Homological_discrete_vector_field::Zp<p, T, false> Type;

    class Is_invertible //   AF: Does not yet exist in Number_types and Algebraic_foundations
      : public CGAL::cpp98::unary_function< Type, bool > {
      public:
        bool operator()( const Type& t) const {
          return t.is_invertible() ;
        }
    };
  };


  // Specialization for p being not a prime number
  template <int p, typename T> class Algebraic_structure_traits< Homological_discrete_vector_field::Zp<p, T, true> >
  : public Algebraic_structure_traits_base< Homological_discrete_vector_field::Zp<p, T, true>, Field_tag >  {
  public:
      typedef Tag_true            Is_exact;
      typedef Tag_false           Is_numerical_sensitive;
      typedef Homological_discrete_vector_field::Zp<p, T, true> Type;

    class Is_invertible
      : public CGAL::cpp98::unary_function< Type, bool > {
      public:
        bool operator()( const Type& t) const {
          return ! t.is_zero() ;
        }
    };
  };

template <int p, typename T, bool IsPrime> class Real_embeddable_traits< Homological_discrete_vector_field::Zp<p, T, IsPrime> >
  : public INTERN_RET::Real_embeddable_traits_base< Homological_discrete_vector_field::Zp<p, T, IsPrime> , CGAL::Tag_true > {
      typedef Homological_discrete_vector_field::Zp<p, T, IsPrime> Type;
  public:

    class Is_positive
      : public CGAL::cpp98::unary_function< Type, bool > {
      public:
        bool operator()( const Type& t) const {
          return ! t.is_zero() ;
        }
    };

    class Is_negative
      : public CGAL::cpp98::unary_function< Type, bool > {
      public:
        bool operator()( const Type&) const {
          return false;
        }
    };

    class Sgn
      : public CGAL::cpp98::unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& t ) const {
          if(t.is_zero()) return ZERO ;
          else return POSITIVE ;
        }
    };

    class Abs
    : public CGAL::cpp98::unary_function< Type, Type > {
      public:
        Type  operator()( const Type& t ) const {
          return t ;
        }
    };
  };

} /* end namespace CGAL */

#endif // CGAL_HDVF_ZP_H
