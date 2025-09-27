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

 The class `Zp` implements the concept `IntegralDomainWithoutDivision` with the field \f$\mathbb Z/p\mathbb Z\f$. This is a "lightweight" implementation which aims at providing fast operations and constructors.

 According to the value of `p`, users can chose the size of the representation used to store values (default size: `int`).

 \warning For \f$\mathbb Z/2\mathbb Z\f$, prefer the class `Z2` which is optimized.

 \cgalModels{IntegralDomainWithoutDivision}

 \tparam p an integer.
 \tparam _TSlot a type used for the inner storage of the values (default: `int`).
 */

template <int p, typename _TSlot = int, bool IsPrime = true>
class Zp {
    _TSlot _i ;
public:

    /** \brief Constructor from a value */
    Zp(_TSlot i=0) : _i( (i>=0)?(i % p):((i % p) + p) ) { }


    // Copy constructor
    Zp(const Zp& a) : _i(a._i) {}


    bool is_zero() const { return _i == 0 ; }

    /** \brief Unary operator+ */
    friend Zp operator+ (const Zp& a)
    {
        return Zp<p, _TSlot>(a) ;
    }

    /** \brief Unary operator-. */
    friend Zp     operator- (const Zp& a)
    {
        return Zp<p, _TSlot>(- a._i) ;
    }

    /** \brief Operator+. */
    friend Zp     operator+ (const Zp& a, const Zp& b)
    {
        return Zp<p, _TSlot, IsPrime>((a._i + b._i)) ;
    }

    /** \brief Operator-. */
    friend Zp     operator- (const Zp& a, const Zp& b)
    {
        return Zp<p, _TSlot, IsPrime>((a._i - b._i)) ;
    }

    /** \brief Operator*. */
    friend Zp     operator* (const Zp& a, const Zp& b)
    {
        return Zp<p, _TSlot, IsPrime>((a._i * b._i)) ;
    }

    /** \brief Operator/. */
    friend Zp     operator/ (const Zp& a, const Zp& b)
    {
        return Zp<p, _TSlot, IsPrime>(a._i / b._i) ;
    }

    /** \brief Operator+=. */
    Zp &     operator+= (const Zp& a)
    {
        _i += a._i ;
        if (_i >= 0)
            _i %= p ;
        else
        {
            _i %= p ;
            _i += p ;
        }
        return *this ;
    }

    /** \brief Operator-=. */
    Zp &     operator-= (const Zp& a)
    {
        _i -= a._i ;
        if (_i >= 0)
            _i %= p ;
        else
        {
            _i %= p ;
            _i += p ;
        }
        return *this ;
    }

    /** \brief Operator*=. */
    Zp &     operator*= (const Zp& a)
    {
        _i *= a._i ;
        if (_i >= 0)
            _i %= p ;
        else
        {
            _i %= p ;
            _i += p ;
        }
        return *this ;
    }

    /** \brief Operator/=. */
    Zp &     operator/= (const Zp& a)
    {
        _i /= a._i ;
        if (_i >= 0)
            _i %= p ;
        else
        {
            _i %= p ;
            _i += p ;
        }
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
        return Zp<p,_TSlot, IsPrime>(abs(int(a._i))) ;
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
};

} /* end namespace Homological_discrete_vector_field */

// Specialization for p being not a prime number
template <int p, typename _TSlot> class Algebraic_structure_traits< Homological_discrete_vector_field::Zp<p, _TSlot, false> >
  : public Algebraic_structure_traits_base< Homological_discrete_vector_field::Zp<p, _TSlot, false>, Integral_domain_without_division_tag >  {
  public:
      typedef Tag_true            Is_exact;
      typedef Tag_false           Is_numerical_sensitive;
      typedef Homological_discrete_vector_field::Zp<p, _TSlot, false> Type;

    class Is_invertible //   AF: Does not yet exist in Number_types and Algebraic_foundations
      : public CGAL::cpp98::unary_function< Type, bool > {
      public:
        bool operator()( const Type& t) const {
          return t.is_invertible() ;
        }
    };
  };


  // Specialization for p being not a prime number
  template <int p, typename _TSlot> class Algebraic_structure_traits< Homological_discrete_vector_field::Zp<p, _TSlot, true> >
  : public Algebraic_structure_traits_base< Homological_discrete_vector_field::Zp<p, _TSlot, true>, Field_tag >  {
  public:
      typedef Tag_true            Is_exact;
      typedef Tag_false           Is_numerical_sensitive;
      typedef Homological_discrete_vector_field::Zp<p, _TSlot, true> Type;

    class Is_invertible
      : public CGAL::cpp98::unary_function< Type, bool > {
      public:
        bool operator()( const Type& t) const {
          return ! t.is_zero() ;
        }
    };
  };

template <int p, typename _TSlot, bool IsPrime> class Real_embeddable_traits< Homological_discrete_vector_field::Zp<p, _TSlot, IsPrime> >
  : public INTERN_RET::Real_embeddable_traits_base< Homological_discrete_vector_field::Zp<p, _TSlot, IsPrime> , CGAL::Tag_true > {
      typedef Homological_discrete_vector_field::Zp<p, _TSlot, IsPrime> Type;
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
