// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: Some comments are original EXACUS comments and aren't adapted. So
//         they may be wrong now.

// TODO: should exponent type be long or Integer ?

#ifndef CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H

#include <CGAL/basic.h>

#if CGAL_USE_LEDA
#include <CGAL/leda_bigfloat.h>
#endif

#if CGAL_USE_CORE
#include <CGAL/CORE_BigFloat.h>
#endif

#if CGAL_USE_MPFR
#include <CGAL/Gmpfr.h>
#endif

#include <CGAL/ipower.h>


namespace CGAL {

namespace internal {

// Don't define default, results in more convinient compiler messages
template< class Type > class Float_traits;
// {
// public:
//   typedef Null_functor    Get_mantissa;
//   typedef Null_functor    Get_exponent;
//   typedef Null_functor    Mul_by_pow_of_2;
// };

#ifdef CGAL_USE_LEDA

// Specialization for leda_bigfloat
template<>
class Float_traits< leda_bigfloat > {
public:
  struct Get_mantissa
    : public CGAL::cpp98::unary_function< leda_bigfloat, leda_integer > {
    leda_integer operator()( const leda_bigfloat& x ) const {
      //std::cout << x.get_significant() << std::endl;
      return x.get_significant();
    }
  };

  struct Get_exponent
    : public CGAL::cpp98::unary_function< leda_bigfloat, long > {
    long operator()( const leda_bigfloat& x ) const {
      return x.get_exponent().to_long();
    }
  };

  struct Mul_by_pow_of_2
    : public CGAL::cpp98::binary_function< leda_bigfloat, long, leda_bigfloat> {
    leda_bigfloat operator()( const leda_bigfloat& a, long e ) const {
      return leda_bigfloat(a.get_significant(), a.get_exponent()+e);
    }
  };
};

#endif

#ifdef CGAL_USE_CORE

// Specialization for CORE::BigFloat
template<>
class Float_traits< CORE::BigFloat > {
public:

  struct Get_mantissa
    : public CGAL::cpp98::unary_function< CORE::BigFloat, CORE::BigInt > {
    CORE::BigInt operator()( const CORE::BigFloat& x ) const {
      return x.m();
    }
  };

  struct Get_exponent
    : public CGAL::cpp98::unary_function< CORE::BigFloat, long > {
    long operator()( const CORE::BigFloat& x ) const {
      return CORE::CHUNK_BIT*x.exp(); // The basis is 2^CORE::CHUNK_BIT
    }
  };

  struct Mul_by_pow_of_2
    : public CGAL::cpp98::binary_function
    < CORE::BigFloat, long , CORE::BigFloat> {
    CORE::BigFloat operator()( const CORE::BigFloat& a, long e ) const {
      return a*CORE::BigFloat::exp2(e);
    }
  };

};
#endif


#if CGAL_USE_MPFR
template<> class Float_traits< Gmpfr > {

  struct Get_mantissa_exponent
    : public CGAL::cpp98::unary_function< Gmpfr, std::pair<Gmpz,long> > {

    std::pair<Gmpz,long> operator()( const Gmpfr& x ) const {
      return x.to_integer_exp();
    }
  };
public:
  struct Get_mantissa
    : public CGAL::cpp98::unary_function< Gmpfr, Gmpz > {
    Gmpz operator()( const Gmpfr& x ) const {
      return Get_mantissa_exponent()(x).first;
    }
  };

  struct Get_exponent
    : public CGAL::cpp98::unary_function< Gmpfr, long > {
    long operator()( const Gmpfr& x ) const {
      return Get_mantissa_exponent()(x).second;
    }
  };

struct Mul_by_pow_of_2
  : public CGAL::cpp98::binary_function< Gmpfr, Gmpz, Gmpfr> {
  Gmpfr operator()( const Gmpfr& a, long e ) const {
    Gmpfr result(0,a.get_precision()); // just to get the prec of a
    if (e >= 0 ){
      mpfr_mul_2si (result.fr(), a.fr(), e, mpfr_get_default_rounding_mode());
      //std::cout << "INPUT   : "<< a <<"+" << e << std::endl;
      //std::cout << "result: "<< result << std::endl;
      //std::cout << "TRUTH   : "<< a * CGAL::ipower(Gmpfr(2),e) << std::endl;
      CGAL_postcondition(a * CGAL::ipower(Gmpfr(2),e) == result);
    }
    else{
      mpfr_div_2si (result.fr(), a.fr(), -e, mpfr_get_default_rounding_mode());
      CGAL_postcondition(a / CGAL::ipower(Gmpfr(2),-e) == result);
    }
    return result;
  }
};
};
#endif
} //namespace internal



} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H
