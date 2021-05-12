// Copyright (c) 2003,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_APOLLONIUS_GRAPH_2_COMPARE_QUADRATIC_H
#define CGAL_APOLLONIUS_GRAPH_2_COMPARE_QUADRATIC_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Apollonius_graph_2/basic.h>
#include <CGAL/functions_on_signs.h>

namespace CGAL {

namespace ApolloniusGraph_2 {

#ifdef COMPARATOR_PROFILER
#include <CGAL/Apollonius_graph_2/comparator_profiler.h>
#endif


//--------------------------------------------------------------------
// help functions for the compulation of various quantites

template < class FT >
inline FT
value_of_D(const FT& a, const FT& b, const FT& c)
{
  return CGAL::square(b) - a * c;
}

template < class FT >
inline FT
value_of_J(const FT& a1, const FT& b1, const FT& a2, const FT& b2)
{
  return (a1 * b2 - a2 * b1);
}

template < class FT >
inline FT
value_of_Jp(const FT& b1, const FT& c1, const FT& b2, const FT& c2)
{
  return (b1 * c2 - b2 * c1);
}

template < class FT >
inline FT
value_of_G(const FT& a1, const FT& c1, const FT& a2, const FT& c2)
{
  return (a1 * c2 - a2 * c1);
}

template < class FT >
inline FT
value_of_K(const FT& a1c2, const FT& a2c1, const FT& b1b2)
{
  return a1c2 + a2c1 - FT(2) * b1b2;
}

template < class FT >
inline FT
value_of_K(const FT& a1, const FT& b1, const FT& c1,
           const FT& a2, const FT& b2, const FT& c2)
{
  return c1 * a2 + a1 * c2 - FT(2) * b1 * b2;
}

template < class FT >
inline FT
value_of_R0(const FT&  J, const FT& Jp,
            const FT& a1, const FT& c1,
            const FT& a2, const FT& c2)
{
  return CGAL::square(a1 * c2 - c1 * a2) - FT(4) * J * Jp;
}

template < class FT >
inline FT
value_of_R0(const FT& D1, const FT& D2, const FT& K)
{
  return CGAL::square(K) - FT(4) * D1 * D2;
}

template < class FT >
inline FT
value_of_P4(const FT& D1, const FT& D2, const FT& K)
{
  return FT(4) * D1 * D2 - CGAL::square(K);
}

template < class FT >
inline FT
value_of_D(const FT& a1, const FT& D1, const FT& a2, const FT& D2)
{
  return D1 * CGAL::square(a2) - D2 * CGAL::square(a1);
}

template < class FT >
inline FT
value_of_P3inf(const FT& a1, const FT& b1, const FT& J, const FT& G)
{
  return FT(2) * b1 * J - a1 * G;
}

template < class FT >
inline FT
value_of_P3inf_bis(const FT& a1, const FT& K, const FT& a2, const FT& D1)
{
  return -(a1 * K + FT(2) * a2 * D1);
}

template < class FT >
inline FT
value_of_P3pinf(const FT& c1, const FT& J, const FT& a1, const FT& Jp)
{
  return c1 * J - a1 * Jp;
}

template < class FT >
inline FT
value_of_L(const FT& a1, const FT& c2, const FT& b1, const FT& b2)
{
  return (a1 * c2 - b1 * b2);
}

template < class FT >
inline FT
value_of_Lp(const FT& a2, const FT& c1, const FT& b1, const FT& b2)
{
  return (a2 * c1 - b1 * b2);
}

template < class FT >
inline FT
value_of_A(const FT& b1, const FT& J, const FT& a1, const FT& L)
{
  return b1 * J + a1 * L;
}

template < class FT >
inline FT
value_of_Q3(const FT& b2, const FT& a2, const FT& J, const FT& L)
{
  return (b2 * J + a2 * L);
}

template < class FT >
inline FT
value_of_Q3(const FT& a2, const FT& b2, const FT& J, const FT& G,
            const FT& K)
{
  return FT(2) * b2 * J - a2 * (G - K);
}

template < class FT >
inline FT
value_of_Q3p(const FT& a1, const FT& b1, const FT& J, const FT& G,
             const FT& K)
{
  return a1 * (G + K) - FT(2) * b1 * J;
}


//--------------------------------------------------------------------

// the trivial method that uses square roots

template < class FT >
inline
Comparison_result
sqrt_compare_l1_l2(const FT& a1, const FT& b1, const FT& c1,
                   const FT& a2, const FT& b2, const FT& c2)
{
  FT D1 = value_of_D(a1, b1, c1);
  FT D2 = value_of_D(a2, b2, c2);

  FT l1 = ( b1 - CGAL::sqrt(D1) ) / a1;
  FT l2 = ( b2 - CGAL::sqrt(D2) ) / a2;

  return CGAL::compare(l1, l2);
}

template < class FT >
inline
Comparison_result
sqrt_compare_l1_r2(const FT& a1, const FT& b1, const FT& c1,
                   const FT& a2, const FT& b2, const FT& c2)
{
  FT D1 = value_of_D(a1, b1, c1);
  FT D2 = value_of_D(a2, b2, c2);

  FT l1 = ( b1 - CGAL::sqrt(D1) ) / a1;
  FT r2 = ( b2 + CGAL::sqrt(D2) ) / a2;

  return CGAL::compare(l1, r2);
}


template < class FT >
inline
Comparison_result
sqrt_compare_r1_l2(const FT& a1, const FT& b1, const FT& c1,
                   const FT& a2, const FT& b2, const FT& c2)
{
  FT D1 = value_of_D(a1, b1, c1);
  FT D2 = value_of_D(a2, b2, c2);

  FT r1 = ( b1 + CGAL::sqrt(D1) ) / a1;
  FT l2 = ( b2 - CGAL::sqrt(D2) ) / a2;

  return CGAL::compare(r1, l2);
}

template < class FT >
inline
Comparison_result
sqrt_compare_r1_r2(const FT& a1, const FT& b1, const FT& c1,
                   const FT& a2, const FT& b2, const FT& c2)
{
  FT D1 = value_of_D(a1, b1, c1);
  FT D2 = value_of_D(a2, b2, c2);

  FT r1 = ( b1 + CGAL::sqrt(D1) ) / a1;
  FT r2 = ( b2 + CGAL::sqrt(D2) ) / a2;

  return CGAL::compare(r1, r2);
}



//--------------------------------------------------------------------


// the DFMT evaluation trees

template < class FT >
inline
Comparison_result
dfmt_compare_l1_l2(const FT& a1, const FT& b1, const FT& c1,
                   const FT& a2, const FT& b2, const FT& c2)
{
  FT J = value_of_J(a1, b1, a2, b2);
  FT K = value_of_K(a1, b1, c1, a2, b2, c2);

  if ( CGAL::is_positive(J) ) {
    if ( CGAL::is_positive(K) )  return SMALLER;  // l1 < l2

    FT D1 = value_of_D(a1, b1, c1);
    FT D2 = value_of_D(a2, b2, c2);

    FT D = value_of_D(a1, D1, a2, D2);

    if ( CGAL::is_positive(D) )  return SMALLER;  // l1 < l2

    FT Jp = value_of_Jp(b1, c1, b2, c2);

    if ( CGAL::is_negative(Jp) )  return LARGER;   // l1 > l2

    FT R0 = value_of_R0(D1, D2, K);

    Sign s_R0 = CGAL::sign(R0);
    if ( s_R0 == NEGATIVE )  return SMALLER;  // l1 < l2
    if ( s_R0 == POSITIVE )  return LARGER;   // l1 > l2
    return EQUAL;
  } else { // J<0
    if ( CGAL::is_positive(K) )  return LARGER;   // l1 > l2

    FT D1 = value_of_D(a1, b1, c1);
    FT D2 = value_of_D(a2, b2, c2);

    FT D = value_of_D(a1, D1, a2, D2);

    if ( CGAL::is_negative(D) )  return LARGER;   // l1 > l2

    FT Jp = value_of_Jp(b1, c1, b2, c2);

    if ( CGAL::is_positive(Jp) )  return SMALLER;  // l1 < l2

    FT R0 = value_of_R0(D1, D2, K);

    Sign s_R0 = CGAL::sign(R0);
    if ( s_R0 == NEGATIVE )  return LARGER;   // l1 > l2
    if ( s_R0 == POSITIVE )  return SMALLER;  // l1 < l2
    return EQUAL;
  }
}

template < class FT >
inline
Comparison_result
dfmt_compare_l1_r2(const FT& a1, const FT& b1, const FT& c1,
                   const FT& a2, const FT& b2, const FT& c2)
{
  FT J = value_of_J(a1, b1, a2, b2);

  if ( CGAL::is_positive(J) ) return SMALLER;   // l1 < r2

  FT K = value_of_K(a1, b1, c1, a2, b2, c2);

  if ( CGAL::is_negative(K) ) return SMALLER;   // l1 < r2

  FT Jp = value_of_Jp(b1, c1, b2, c2);

  if ( CGAL::is_positive(Jp) ) return LARGER;  // l1 > r2

  FT R0 = value_of_R0(J, Jp, a1, c1, a2, c2);

  Sign s_R0 = CGAL::sign(R0);
  if ( s_R0 == NEGATIVE ) return SMALLER;   // l1 < r2
  if ( s_R0 == POSITIVE ) return LARGER;    // l1 > r2
  return EQUAL;
}



template < class FT >
inline
Comparison_result
dfmt_compare_r1_l2(const FT& a1, const FT& b1, const FT& c1,
                   const FT& a2, const FT& b2, const FT& c2)
{
  FT J = value_of_J(a1, b1, a2, b2);

  if ( CGAL::is_negative(J) ) return LARGER;   // r1 > l2

  FT K = value_of_K(a1, b1, c1, a2, b2, c2);

  if ( CGAL::is_negative(K) ) return LARGER;   // r1 > l2

  FT Jp = value_of_Jp(b1, c1, b2, c2);

  if ( CGAL::is_negative(Jp) ) return SMALLER;  // r1 < l2

  FT R0 = value_of_R0(J, Jp, a1, c1, a2, c2);

  Sign s_R0 = CGAL::sign(R0);
  if ( s_R0 == NEGATIVE ) return LARGER;   // r1 > l2
  if ( s_R0 == POSITIVE ) return SMALLER;  // r1 < l2
  return EQUAL;
}

template < class FT >
inline
Comparison_result
dfmt_compare_r1_r2(const FT& a1, const FT& b1, const FT& c1,
                   const FT& a2, const FT& b2, const FT& c2)
{
#ifdef COMPARATOR_PROFILER
  comparator_profiler::counter_rr++;
#endif

  FT J = value_of_J(a1, b1, a2, b2);
  FT K = value_of_K(a1, b1, c1, a2, b2, c2);

  if ( CGAL::is_positive(J) ){
    if ( CGAL::is_positive(K) )  return SMALLER;   // r1 < r2   1,2

#ifdef COMPARATOR_PROFILER
    comparator_profiler::counter_rr_e++;
#endif

    FT D1 = value_of_D(a1, b1, c1);
    FT D2 = value_of_D(a2, b2, c2);

    FT D = value_of_D(a1, D1, a2, D2);

    if ( CGAL::is_negative(D) )  return SMALLER;   // r1 < r2   2,3b

#ifdef COMPARATOR_PROFILER
    comparator_profiler::counter_rr_r0++;
#endif

    FT Jp = value_of_Jp(b1, c1, b2, c2);

    if ( CGAL::is_negative(Jp) )  return LARGER;    // r1 > r2   3a

    FT R0 = value_of_R0(D1, D2, K);

    Sign s_R0 = CGAL::sign(R0);
    if ( s_R0 == NEGATIVE )  return SMALLER;   // r1 < r2   2
    if ( s_R0 == POSITIVE )  return LARGER;    // r1 > r2   3a
    return EQUAL;
  } else { // J<0
    if ( CGAL::is_positive(K) )  return LARGER;   // r1 > r2    4,5

#ifdef COMPARATOR_PROFILER
    comparator_profiler::counter_rr_e++;
#endif

    FT D1 = value_of_D(a1, b1, c1);
    FT D2 = value_of_D(a2, b2, c2);

    FT D = value_of_D(a1, D1, a2, D2);

    if ( CGAL::is_positive(D) )  return LARGER;   // r1 > r2    3a,4

#ifdef COMPARATOR_PROFILER
    comparator_profiler::counter_rr_r0++;
#endif

    FT Jp = value_of_Jp(b1, c1, b2, c2);

    if ( CGAL::is_positive(Jp) )  return SMALLER;   // r1 < r2     3b

    FT R0 = value_of_P4(D1, D2, K);

    Sign s_R0 = CGAL::sign(R0);
    if ( s_R0 == NEGATIVE )  return LARGER;   // r1 > r2     4
    if ( s_R0 == POSITIVE )  return SMALLER;  // r1 < r2     3b
    return EQUAL;
  }
}

//--------------------------------------------------------------------

// the KE evaluation trees

template < class FT >
inline
Comparison_result
ke_compare_l1_l2(const FT& a1, const FT& b1, const FT& c1,
                 const FT& a2, const FT& b2, const FT& c2)
{
  FT J = value_of_J(a1, b1, a2, b2);
  Sign s_J = CGAL::sign(J);

  if ( s_J == ZERO ) {
    Sign s_G = CGAL::sign( value_of_G(a1, c1, a2, c2) );
    if ( s_G == POSITIVE ) { return SMALLER; }
    if ( s_G == NEGATIVE ) { return LARGER; }
    return EQUAL;
  }

  FT a1c2 = a1 * c2;
  FT a2c1 = a2 * c1;
  FT K = value_of_K<FT>(a1c2, a2c1, b1 * b2);
  Sign s_K = CGAL::sign(K);

  if ( s_J == POSITIVE ) {
    if ( s_K == POSITIVE ) { return SMALLER; }

    if ( s_K == ZERO ) {
      FT D1 = value_of_D(a1, b1, c1);
      if ( CGAL::is_zero(D1) ) { return EQUAL; }

      return SMALLER;
    }

    FT G = a1c2 - a2c1;
    FT P3inf = value_of_P3inf(a1, b1, J, G);

    if ( !(CGAL::is_positive(P3inf)) ) { return SMALLER; }

    FT Jp = value_of_Jp(b1, c1, b2, c2);

    if ( CGAL::is_negative(Jp) )  { return LARGER; }

    FT P4 = value_of_P4(J, Jp, G);

    Sign s_P4 = CGAL::sign(P4);
    if ( s_P4 == POSITIVE )  { return SMALLER; }
    if ( s_P4 == NEGATIVE )  { return LARGER; }
    return EQUAL;
  }

  // J < 0
  if ( s_K == POSITIVE )  { return LARGER; }

  if ( s_K == ZERO ) {
    FT D2 = value_of_D(a2, b2, c2);
    if ( CGAL::is_zero(D2) ) { return EQUAL; }

    return LARGER;
  }
  FT G = a1c2 - a2c1;
  FT P3inf = value_of_P3inf(a1, b1, J, G);

  if ( !(CGAL::is_negative(P3inf)) )  { return LARGER; }

  FT Jp = value_of_Jp(b1, c1, b2, c2);

  if ( CGAL::is_positive(Jp) )  { return SMALLER; }

  FT P4 = value_of_P4(J, Jp, G);

  Sign s_P4 = CGAL::sign(P4);
  if ( s_P4 == POSITIVE )  { return LARGER; }
  if ( s_P4 == NEGATIVE )  { return SMALLER; }
  return EQUAL;
}

template < class FT >
inline
Comparison_result
ke_compare_l1_r2(const FT& a1, const FT& b1, const FT& c1,
                 const FT& a2, const FT& b2, const FT& c2)
{
  FT J = value_of_J(a1, b1, a2, b2);
  Sign s_J = CGAL::sign(J);

#if 0
  if ( s_J == ZERO ) {
    FT D1 = value_of_D(a1, b1, c1);
    if ( CGAL::is_positive(D1) ) { return SMALLER; }

    FT D2 = value_of_D(a2, b2, c2);
    if ( CGAL::is_positive(D2) ) { return SMALLER; }

    return EQUAL;
  }
#endif

  if ( s_J == POSITIVE ) { return SMALLER; }

  FT a1c2 = a1 * c2;
  FT a2c1 = a2 * c1;
  FT K = value_of_K<FT>(a1c2, a2c1, b1 * b2);
  Sign s_K = CGAL::sign(K);

  if ( s_K == NEGATIVE ) { return SMALLER; }
#if 0
  if ( s_K == ZERO ) {
    FT D1 = value_of_D(a1, b1, c1);
    if ( CGAL::is_zero(D1) ) { return EQUAL; }

    FT D2 = value_of_D(a2, b2, c2);
    if ( CGAL::is_zero(D2) ) { return EQUAL; }

    return SMALLER;
  }
#endif

  FT Jp = value_of_Jp(b1, c1, b2, c2);

  if ( CGAL::is_positive(Jp) ) { return LARGER; }

  FT P4 = value_of_P4<FT>(J, Jp, a1c2 - a2c1);

  Sign s_P4 = CGAL::sign(P4);
  if ( s_P4 == POSITIVE ) { return SMALLER; }
  if ( s_P4 == NEGATIVE ) { return LARGER; }
  return EQUAL;
}

template < class FT >
inline
Comparison_result
ke_compare_r1_l2(const FT& a1, const FT& b1, const FT& c1,
                 const FT& a2, const FT& b2, const FT& c2)
{
  FT J = value_of_J(a1, b1, a2, b2);
  Sign s_J = CGAL::sign(J);

#if 0
  if ( s_J == ZERO ) {
    FT D1 = value_of_D(a1, b1, c1);
    if ( CGAL::is_positive(D1) ) { return LARGER; }

    FT D2 = value_of_D(a2, b2, c2);
    if ( CGAL::is_positive(D2) ) { return LARGER; }

    return EQUAL;
  }
#endif

  if ( s_J == NEGATIVE ) { return LARGER; }

  FT a1c2 = a1 * c2;
  FT a2c1 = a2 * c1;
  FT K = value_of_K<FT>(a1c2, a2c1, b1 * b2);
  Sign s_K = CGAL::sign(K);

  if ( s_K == NEGATIVE ) { return LARGER; }

#if 0
  if ( s_K == ZERO ) {
    FT D1 = value_of_D(a1, b1, c1);
    if ( CGAL::is_zero(D1) ) { return EQUAL; }

    FT D2 = value_of_D(a2, b2, c2);
    if ( CGAL::is_zero(D2) ) { return EQUAL; }

    return LARGER;
  }
#endif

  FT Jp = value_of_Jp(b1, c1, b2, c2);

  if ( CGAL::is_negative(Jp) ) { return SMALLER; }

  FT P4 = value_of_P4<FT>(J, Jp, a1c2 - a2c1);

  Sign s_P4 = CGAL::sign(P4);
  if ( s_P4 == POSITIVE ) { return LARGER; }
  if ( s_P4 == NEGATIVE ) { return SMALLER; }
  return EQUAL;
}

template < class FT >
inline
Comparison_result
ke_compare_r1_r2(const FT& a1, const FT& b1, const FT& c1,
                 const FT& a2, const FT& b2, const FT& c2)
{
#ifdef COMPARATOR_PROFILER
  comparator_profiler::counter_rr++;
#endif

  FT J = value_of_J(a1, b1, a2, b2);
  Sign s_J = CGAL::sign(J);

  FT a1c2 = a1 * c2;
  FT a2c1 = a2 * c1;
  FT K = value_of_K<FT>(a1c2, a2c1, b1 * b2);
  Sign s_K = CGAL::sign(K);

  if ( s_J == POSITIVE ) {
    if ( s_K == POSITIVE )  { return SMALLER; }
    else if ( s_K == NEGATIVE ) {

#ifdef COMPARATOR_PROFILER
      comparator_profiler::counter_rr_p3inf++;
#endif

      FT G = a1c2 - a2c1;
      FT P3inf = value_of_P3inf(a1, b1, J, G);

      if ( !(CGAL::is_negative(P3inf)) ) { return SMALLER; }

#ifdef COMPARATOR_PROFILER
      comparator_profiler::counter_rr_p4++;
#endif

      FT Jp = value_of_Jp(b1, c1, b2, c2);

      if ( CGAL::is_negative(Jp) )  { return LARGER; }

      FT P4 = value_of_P4(J, Jp, G);

      Sign s_P4 = CGAL::sign(P4);
      if ( s_P4 == POSITIVE )  { return SMALLER; }
      if ( s_P4 == NEGATIVE )  { return LARGER; }
      return EQUAL;
    } else { // K = 0
      FT D2 = value_of_D(a2, b2, c2);
      if ( CGAL::is_zero(D2) ) { return EQUAL; }

      return SMALLER;
    }
  } else if ( s_J == NEGATIVE ) { // J < 0
    if ( s_K == POSITIVE )  { return LARGER; }
    else if ( s_K == NEGATIVE ) {

#ifdef COMPARATOR_PROFILER
      comparator_profiler::counter_rr_p3inf++;
#endif

      FT G = a1c2 - a2c1;
      FT P3inf = value_of_P3inf(a1, b1, J, G);

      if ( !(CGAL::is_positive(P3inf)) ) { return LARGER; }

#ifdef COMPARATOR_PROFILER
      comparator_profiler::counter_rr_p4++;
#endif

      FT Jp = value_of_Jp(b1, c1, b2, c2);

      if ( CGAL::is_positive(Jp) )   { return SMALLER; }

      FT P4 = value_of_P4(J, Jp, G);

      Sign s_P4 = CGAL::sign(P4);
      if ( s_P4 == POSITIVE )  { return LARGER; }
      if ( s_P4 == NEGATIVE )  { return SMALLER; }
      return EQUAL;
    } else { // K = 0
      FT D1 = value_of_D(a1, b1, c1);
      if ( CGAL::is_zero(D1) ) { return EQUAL; }

      return LARGER;
    }
  }

  // J = 0
  Sign s_G = CGAL::sign( value_of_G(a1, c1, a2, c2) );
  if ( s_G == NEGATIVE ) { return SMALLER; }
  if ( s_G == POSITIVE ) { return LARGER; }
  return EQUAL;
}


//--------------------------------------------------------------------

//--------------------------------------------------------------------
//--------------------------------------------------------------------
// the following functions do the filtering for the r1-r2 tree without
// using C++ exceptions

template < class CT, class ET >
inline
Comparison_result
sqrt_compare_r1_r2_filtered(const CT& a1, const CT& b1, const CT& c1,
                            const CT& a2, const CT& b2, const CT& c2)
{
  typedef Interval_nt<false> IT;

  FPU_CW_t backup = FPU_get_cw();
  FPU_set_cw(CGAL_FE_UPWARD);

  IT a1i(a1), b1i(b1), c1i(c1);
  IT a2i(a2), b2i(b2), c2i(c2);

  IT D1 = value_of_D(a1i, b1i, c1i);
  IT D2 = value_of_D(a2i, b2i, c2i);

  IT r1 = ( b1i + CGAL::sqrt(D1) ) / a1i;
  IT r2 = ( b2i + CGAL::sqrt(D2) ) / a2i;

  FPU_set_cw(backup);

  if ( r1.sup() < r2.inf() ) { return SMALLER; }
  if ( r1.inf() > r2.sup() ) { return LARGER; }

  return sqrt_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                            ET(a2), ET(b2), ET(c2));
}

//--------------------------------------------------------------------

template < class CT, class ET >
inline
Comparison_result
dfmt_compare_r1_r2_filtered(const CT& a1, const CT& b1, const CT& c1,
                            const CT& a2, const CT& b2, const CT& c2)
{
  typedef Interval_nt<false> IT;

  FPU_CW_t backup = FPU_get_cw();
  FPU_set_cw(CGAL_FE_UPWARD);

  IT a1i(a1), b1i(b1), c1i(c1);
  IT a2i(a2), b2i(b2), c2i(c2);

  IT J = value_of_J(a1i, b1i, a2i, b2i);
  IT K = value_of_K(a1i, b1i, c1i, a2i, b2i, c2i);

  if ( J.inf() > 0 ){
    if ( K.inf() > 0 ) {
      FPU_set_cw(backup);
      return SMALLER;
    }
    if ( K.sup() > 0 ) {
      FPU_set_cw(backup);
      return dfmt_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                                ET(a2), ET(b2), ET(c2));
    }

    IT D1 = value_of_D(a1i, b1i, c1i);
    IT D2 = value_of_D(a2i, b2i, c2i);

    IT D = value_of_D(a1i, D1, a2i, D2);

    if ( D.sup() < 0 ) {
      FPU_set_cw(backup);
      return SMALLER;
    }
    if ( D.inf() < 0 ) {
      FPU_set_cw(backup);
      return dfmt_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                                ET(a2), ET(b2), ET(c2));
    }

    IT Jp = value_of_Jp(b1i, c1i, b2i, c2i);

    if ( Jp.sup() < 0 ) {
      FPU_set_cw(backup);
      return LARGER;
    }
    if ( Jp.inf() < 0 ) {
      FPU_set_cw(backup);
      return dfmt_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                                ET(a2), ET(b2), ET(c2));
    }


    IT R0 = value_of_R0(D1, D2, K);
    FPU_set_cw(backup);

    if ( R0.sup() < 0 ) { return SMALLER; }
    if ( R0.inf() > 0 ) { return LARGER; }
    return dfmt_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                              ET(a2), ET(b2), ET(c2));
  } else if ( J.sup() > 0 ) {
    FPU_set_cw(backup);
    return dfmt_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                              ET(a2), ET(b2), ET(c2));
  } else { // J < 0
    if ( K.inf() > 0 ) {
      FPU_set_cw(backup);
      return LARGER;
    }
    if ( K.sup() > 0 ) {
      FPU_set_cw(backup);
      return dfmt_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                                ET(a2), ET(b2), ET(c2));
    }

    IT D1 = value_of_D(a1i, b1i, c1i);
    IT D2 = value_of_D(a2i, b2i, c2i);

    IT D = value_of_D(a1i, D1, a2i, D2);

    if ( D.inf() > 0 ) {
      FPU_set_cw(backup);
      return LARGER;
    }
    if ( D.sup() > 0 ) {
      FPU_set_cw(backup);
      return dfmt_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                                ET(a2), ET(b2), ET(c2));
    }


    IT Jp = value_of_Jp(b1i, c1i, b2i, c2i);

    if ( Jp.inf() > 0 ) {
       FPU_set_cw(backup);
      return SMALLER;
    }
    if ( Jp.sup() > 0 ) {
      FPU_set_cw(backup);
      return dfmt_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                                ET(a2), ET(b2), ET(c2));
    }

    IT R0 = value_of_P4(D1, D2, K);
    FPU_set_cw(backup);

    if ( R0.sup() < 0 ) { return LARGER; }
    if ( R0.inf() > 0 ) { return SMALLER; }

    return dfmt_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                              ET(a2), ET(b2), ET(c2));
  }
}

//--------------------------------------------------------------------

template < class CT, class ET >
inline
Comparison_result
ke_compare_r1_r2_filtered(const CT& a1, const CT& b1, const CT& c1,
                          const CT& a2, const CT& b2, const CT& c2)
{
  typedef Interval_nt<false> IT;

  FPU_CW_t backup = FPU_get_cw();
  FPU_set_cw(CGAL_FE_UPWARD);

  IT a1i(a1), b1i(b1), c1i(c1);
  IT a2i(a2), b2i(b2), c2i(c2);

  IT J = value_of_J(a1i, b1i, a2i, b2i);

  IT a1c2 = a1i * c2i;
  IT a2c1 = a2i * c1i;

  IT K = value_of_K(a1c2, a2c1, b1i * b2i);

  if ( J.inf() > 0 ) {
    if ( K.inf() > 0 )  {
      FPU_set_cw(backup);
      return SMALLER;
    } else if ( K.sup() < 0 ) {
      IT G = a1c2 - a2c1;
      IT P3inf = value_of_P3inf(a1i, b1i, J, G);

      if ( P3inf.inf() >= 0 ) {
        FPU_set_cw(backup);
        return SMALLER;
      }
      if ( P3inf.sup() > 0 ) {
        FPU_set_cw(backup);
        return ke_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                                ET(a2), ET(b2), ET(c2));
      }

      IT Jp = value_of_Jp(b1i, c1i, b2i, c2i);

      if ( Jp.sup() < 0 ) {
        FPU_set_cw(backup);
        return LARGER;
      }
      if ( Jp.inf() < 0 ) {
        FPU_set_cw(backup);
        return ke_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                                ET(a2), ET(b2), ET(c2));
      }

      IT P4 = value_of_P4(J, Jp, G);
      FPU_set_cw(backup);

      if ( P4.inf() > 0 ) { return SMALLER; }
      if ( P4.sup() < 0 ) { return LARGER; }

      return ke_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                              ET(a2), ET(b2), ET(c2));
    } else {
      FPU_set_cw(backup);
      return ke_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                              ET(a2), ET(b2), ET(c2));
    }
  } else if ( J.sup() < 0 ) { // J < 0
    if ( K.inf() > 0 ) {
      FPU_set_cw(backup);
      return LARGER;
    } else if ( K.sup() < 0 ) {
      IT G = a1c2 - a2c1;
      IT P3inf = value_of_P3inf(a1i, b1i, J, G);

      if ( P3inf.sup() <= 0 ) {
        FPU_set_cw(backup);
        return LARGER;
      }
      if ( P3inf.inf() < 0 ) {
        FPU_set_cw(backup);
        return ke_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                                ET(a2), ET(b2), ET(c2));
      }

      IT Jp = value_of_Jp(b1i, c1i, b2i, c2i);

      if ( Jp.inf() > 0 ) {
        FPU_set_cw(backup);
        return SMALLER;
      }
      if ( Jp.sup() > 0 ) {
        FPU_set_cw(backup);
        return ke_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                                ET(a2), ET(b2), ET(c2));
      }

      IT P4 = value_of_P4(J, Jp, G);
      FPU_set_cw(backup);

      if ( P4.inf() > 0 ) { return LARGER; }
      if ( P4.sup() < 0 ) { return SMALLER; }

      return ke_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                              ET(a2), ET(b2), ET(c2));
    } else {
      FPU_set_cw(backup);
      return ke_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                              ET(a2), ET(b2), ET(c2));
    }
  } else { // J = ?
    FPU_set_cw(backup);
    return ke_compare_r1_r2(ET(a1), ET(b1), ET(c1),
                            ET(a2), ET(b2), ET(c2));
  }
}





//--------------------------------------------------------------------
//--------------------------------------------------------------------

} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_COMPARE_QUADRATIC_H
