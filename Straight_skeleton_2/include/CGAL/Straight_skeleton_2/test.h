// Copyright (c) 2007 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>

#ifndef CGAL_STRAIGHT_SKELETON_TEST_H
#define CGAL_STRAIGHT_SKELETON_TEST_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/number_utils.h>

#include <algorithm>
#include <iostream>

//
// INTRINSIC UNIT TESTING macros.
//
// These are OFF by default and should only be turn on by the testsuite.
//

//
// These are assertions that may not be true unless certain external conditions are met, such as the kernel
// beigng exact or the input being specifically designed to be non-degenerate, avoid significant round-off, etc.
//
#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_INTRINSIC_TESTING
#  define CGAL_stskel_intrinsic_test_assertion(EX)          ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define CGAL_stskel_intrinsic_test_assertion_msg(EX,MSG)  ((EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_stskel_intrinsic_test_assertion_code(CODE)   CODE
#  define CGAL_stskel_intrinsic_test_trace(m)               std::cerr << m << std::endl
#  define CGAL_stskel_intrinsic_test_trace_if(EX,m)         if ( (EX) ) { std::cerr << m << std::endl ; }
#else
#  define CGAL_stskel_intrinsic_test_assertion(EX)         (static_cast<void>(0))
#  define CGAL_stskel_intrinsic_test_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_stskel_intrinsic_test_assertion_code(CODE)
#  define CGAL_stskel_intrinsic_test_trace(m)
#  define CGAL_stskel_intrinsic_test_trace_if(EX,m)
#endif

namespace CGAL {

namespace CGAL_SS_i {

//
// The following tests are used by the testsuite only.
//
// Their purpose is to detect clearly wrong results without resorting to exact constructions.
//
// These are negative tests only. For instance, they don't test whether a number is zero (since it can be
// near zero but not exactly due to roundoff); rather, they test whether a number is clearly not zero,
// which is a test that can be done robustelly if a pesimistic upper bound on the error is know.
//
// The test are overloaded on number types so if exact constructions are used, the tests are exact.

inline bool is_possibly_inexact_distance_clearly_not_zero ( double n, double eps )
{
  return std::abs( CGAL_NTS to_double(n) ) > eps ;
}

#ifdef CGAL_CORE_EXPR_H
inline bool is_possibly_inexact_distance_clearly_not_zero ( CORE::Expr const& n )
{
  return ! CGAL_NTS is_zero(n);
}
#endif

#ifdef CGAL_LEDA_REAL_H
inline bool is_possibly_inexact_distance_clearly_not_zero ( leda_real const& n )
{
  return ! CGAL_NTS is_zero(n);
}
#endif

#ifdef CGAL_GMPQ_H
inline bool is_possibly_inexact_distance_clearly_not_zero ( Gmpq const& n )
{
  return is_possibly_inexact_distance_clearly_not_zero( to_double(n), 1e-8 ) ;
}
#endif

#ifdef CGAL_MP_FLOAT_H
inline bool is_possibly_inexact_distance_clearly_not_zero ( MP_Float const& n )
{
  return is_possibly_inexact_distance_clearly_not_zero( to_double(n), 1e-8 ) ;
}

inline bool is_possibly_inexact_distance_clearly_not_zero ( Quotient<MP_Float> const& n )
{
  return is_possibly_inexact_distance_clearly_not_zero( to_double(n), 1e-8 ) ;
}
#endif

#if defined(CGAL_LAZY_EXACT_NT_H)
template<class NT>
inline bool is_possibly_inexact_distance_clearly_not_zero ( Lazy_exact_nt<NT> const& n )
{
  return is_possibly_inexact_distance_clearly_not_zero( to_double(n), 1e-8 ) ;
}
#endif


inline bool is_possibly_inexact_distance_clearly_not_zero ( double n )
{
  return std::abs( CGAL_NTS to_double(n) ) > 1e-5 ;
}

inline bool is_possibly_inexact_distance_clearly_not_zero ( Interval_nt_advanced const& n )
{
  return is_possibly_inexact_distance_clearly_not_zero(to_double(n));
}





template<class NT>
inline bool is_possibly_inexact_distance_clearly_not_equal_to( NT const& n, NT const& m )
{
  return is_possibly_inexact_distance_clearly_not_zero(n-m);
}

template<class NT>
inline bool is_possibly_inexact_time_clearly_not_zero( NT const& n )
{
  return is_possibly_inexact_distance_clearly_not_zero(n);
}

template<class NT>
inline bool is_possibly_inexact_time_clearly_not_equal_to( NT const& n, NT const& m )
{
  return is_possibly_inexact_distance_clearly_not_zero(n-m);
}

template<class FT, class Bisector>
inline bool is_time_clearly_not_within_possibly_inexact_bisector_time_interval( FT const& aT , Bisector const& aBisector )
{
  FT lSrcT = aBisector->opposite()->vertex()->time() ;
  FT lTgtT = aBisector->vertex()->time() ;
  FT lLoT  = (std::min)(lSrcT,lTgtT);
  FT lHiT  = (std::max)(lSrcT,lTgtT);

  return    ( aT < lLoT || aT > lHiT )
         && is_possibly_inexact_time_clearly_not_equal_to(aT,lLoT)
         && is_possibly_inexact_time_clearly_not_equal_to(aT,lHiT) ;
}

template<class FT, class Bisector>
inline bool is_time_clearly_within_possibly_inexact_bisector_time_interval( FT const& aT , Bisector const& aBisector )
{
  FT lSrcT = aBisector->opposite()->vertex()->time() ;
  FT lTgtT = aBisector->vertex()->time() ;
  FT lLoT  = (std::min)(lSrcT,lTgtT);
  FT lHiT  = (std::max)(lSrcT,lTgtT);

  return    ( lLoT < aT && aT < lHiT )
         && is_possibly_inexact_time_clearly_not_equal_to(aT,lLoT)
         && is_possibly_inexact_time_clearly_not_equal_to(aT,lHiT) ;
}


} // namespace CGAL_SS_i

} // end namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_TEST_H //
// EOF //


