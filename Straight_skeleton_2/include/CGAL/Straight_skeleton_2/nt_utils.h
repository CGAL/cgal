// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>

#ifndef CGAL_STRAIGHT_SKELETON_NT_UTILS_H
#define CGAL_STRAIGHT_SKELETON_NT_UTILS_H 1



#ifdef CGAL_USE_CORE
#  include <CGAL/CORE_BigFloat.h>
#endif

namespace CGAL{

namespace CGAL_SS_i {

#ifdef CGAL_USE_CORE

#ifdef CGAL_USE_GMPXX
inline CORE::BigFloat to_BigFloat( mpq_class const& n )
{
  return CORE::BigFloat( CORE::BigRat(n.get_mpq_t()) );
}
#else
inline CORE::BigFloat to_BigFloat( CGAL::Gmpq const& n )
{
  return CORE::BigFloat( CORE::BigRat(n.mpq()) );
}
#endif

template<class NT>
inline CORE::BigFloat to_BigFloat( NT const& n )
{
  return to_BigFloat( CGAL::exact(n) );
}

template<>
inline CORE::BigFloat to_BigFloat<MP_Float>( MP_Float const& b )
{
  if (b.is_zero())
    return CORE::BigFloat::getZero();

  typedef MP_Float::exponent_type exponent_type;

  const int                    log_limb         = 8 * sizeof(MP_Float::limb);
  const MP_Float::V::size_type limbs_per_double = 2 + 53/log_limb;

  exponent_type exp = b.max_exp();
  int steps = static_cast<int>((std::min)(limbs_per_double, b.v.size()));

  CORE::BigFloat d_exp_1 = CORE::BigFloat::exp2(-log_limb);

  CORE::BigFloat d_exp   = CORE::BigFloat::getOne() ;

  CORE::BigFloat d       = CORE::BigFloat::getZero();

  for ( exponent_type i = exp - 1; i > exp - 1 - steps; i--)
  {
    d_exp *= d_exp_1;
    d += d_exp * CORE::BigFloat(b.of_exp(i));
  }

  return d * CORE::BigFloat::exp2( static_cast<int>(exp * log_limb) );
}


#endif

template<class NT>
inline NT inexact_sqrt_implementation( NT const& n, CGAL::Null_functor /*no_sqrt*/ )
{

#ifdef CGAL_USE_CORE

  CORE::BigFloat nn = to_BigFloat(n) ;
  CORE::BigFloat s  = CORE::sqrt(nn);
  return NT(s.doubleValue());

#else

  double nn = CGAL::to_double(n) ;

  if ( !CGAL_NTS is_valid(nn) || ! CGAL_NTS is_finite(nn) )
    nn = std::numeric_limits<double>::max BOOST_PREVENT_MACRO_SUBSTITUTION () ;

  CGAL_precondition(nn > 0);

  double s = CGAL_NTS sqrt(nn);

  return NT(s);

#endif
}

template<class NT, class Sqrt>
inline NT inexact_sqrt_implementation( NT const& n, Sqrt sqrt_f )
{
  return sqrt_f(n);
}

template<class NT>
inline NT inexact_sqrt( NT const& n )
{
  typedef CGAL::Algebraic_structure_traits<NT> AST;
  typedef typename AST::Sqrt Sqrt;
  return inexact_sqrt_implementation(n,Sqrt());
}

inline Quotient<MP_Float> inexact_sqrt( Quotient<MP_Float> const& q )
{
  CGAL_precondition(q > 0);
  return Quotient<MP_Float>(CGAL_SS_i::inexact_sqrt(q.numerator()*q.denominator()), q.denominator() );
}

} } // end of CGAL::CGAL_SS_i namespace

#endif // CGAL_STRAIGHT_SKELETON_NT_UTILS_H
