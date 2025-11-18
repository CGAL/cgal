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

#ifndef CGAL_HDVF_Z2_H
#define CGAL_HDVF_Z2_H

#include <CGAL/license/HDVF.h>

#include <CGAL/number_type_basic.h>
#include <CGAL/number_utils.h>

#include <vector>
#include <iostream>

namespace CGAL {
namespace Homological_discrete_vector_field {

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Z2` implements the concept `IntegralDomainWithoutDivision` with the field \f$\mathbb Z/2\mathbb Z\f$. This implementation is optimized to use bitwise operations and should be prefered to `Zp<2>`.

 \cgalModels{IntegralDomainWithoutDivision}
 */

class Z2 {
    char _i ;
public:

    /** \brief Constructor from a value */
    Z2(char i=0) : _i(i ? 1 : 0 ) {}

    // Copy constructor
    Z2(const Z2& a) : _i(a._i) {}
    
    /** \brief Returns 2 (ie. p such that `Z2`=\f$\mathbbZ/2\mathbb Z\f$). */
    static char operator() () { return char(2); }

    bool is_zero() const { return _i == 0 ; }

    /** \brief Unary operator+ */
    friend Z2 operator+ (const Z2& a)
    {
        return a ;
    }

    /** \brief Unary operator-. */
    friend Z2     operator- (const Z2& a)
    {
        return a ;
    }

    /** \brief Operator+. */
    friend Z2     operator+ (const Z2& a, const Z2& b)
    {
        return Z2((a._i^b._i)) ;
    }

    /** \brief Operator-. */
    friend Z2     operator- (const Z2& a, const Z2& b)
    {
        return Z2(a._i^b._i) ;
    }

    /** \brief Operator*. */
    friend Z2     operator* (const Z2& a, const Z2& b)
    {
        return Z2(a._i & b._i) ;
    }

    /** \brief Operator/. */
    friend Z2     operator/ (const Z2& a, const Z2& b)
    {
        return Z2(a._i / b._i) ;
    }

    /** \brief Operator+=. */
    Z2 &     operator+= (const Z2& a)
    {
        _i ^= a._i ;
        return *this ;
    }

    /** \brief Operator-=. */
    Z2 &     operator-= (const Z2& a)
    {
        _i ^= a._i ;
        return *this ;
    }

    /** \brief Operator*=. */
    Z2 &     operator*= (const Z2& a)
    {
        _i &= a._i ;
        return *this ;
    }

    /** \brief Operator/=. */
    Z2 &     operator/= (const Z2& a)
    {
        _i /= a._i ;
        return *this ;
    }

    /** \brief Operator==. */
    friend bool     operator== (const Z2& a, const Z2& b)
    {
        return (a._i == b._i) ;
    }

    /** \brief Operator!=. */
    friend bool     operator!= (const Z2& a, const Z2& b)
    {
        return (a._i != b._i);
    }
#if 0
    /** \brief Absolute value. */
    friend  Z2  abs(const Z2& a)
    {
        return a ;
    }
#endif
    /** \brief Operator<<. */
    friend std::ostream& operator<<(std::ostream& out, const Z2& a)
    {
        return (out << (a._i ? 1 : 0)) ;
    }

    /** \brief Operator>>. */
    friend std::istream& operator>>(std::istream& in, Z2& a)
    {
        int tmp ;
        in >> tmp ;
        a = Z2(tmp) ;
        return (in) ;
    }
};



} /* end namespace Homological_discrete_vector_field */

template <> class Algebraic_structure_traits< Homological_discrete_vector_field::Z2 >
  : public Algebraic_structure_traits_base< Homological_discrete_vector_field::Z2, Integral_domain_without_division_tag >  {
  public:
    typedef Tag_true            Is_exact;
    typedef Tag_false           Is_numerical_sensitive;

    class Is_square
      : public CGAL::cpp98::unary_function< Type, bool > {
      public:
        bool operator()( const Type& ) const {
          return true;
        }
    };
  };

  template <> class Real_embeddable_traits< Homological_discrete_vector_field::Z2 >
  : public INTERN_RET::Real_embeddable_traits_base< Homological_discrete_vector_field::Z2 , CGAL::Tag_true > {
  public:

    class Is_positive
      : public CGAL::cpp98::unary_function< Type, bool > {
      public:
        bool operator()( const Type& x_) const {
          return !x_.is_zero() ;
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
        ::CGAL::Sign operator()( const Type& x ) const {
          if(x.is_zero()) return ZERO ;
          else return POSITIVE ;
        }
    };

    class Abs
    : public CGAL::cpp98::unary_function< Type, Type > {
      public:
        Type  operator()( const Type& x ) const {
          return x ;
        }
    };
  };

} /* end namespace CGAL */

#endif //CGAL_HDVF_Z2_H
