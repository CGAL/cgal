// Copyright (c) 2005,2006  INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Sylvain Pion, Monique Teillaud, Athanasios Kakargias
//                 Olivier Devillers

#ifndef CGAL_NUMBER_TypeS_ROOT_OF_COMPARISON_FUNCTIONS_22_H
#define CGAL_NUMBER_TypeS_ROOT_OF_COMPARISON_FUNCTIONS_22_H

#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>

namespace CGAL {
namespace internal {

// Maybe we can trash this
/*1 1*/template <class FT>
/*1 1*/Comparison_result
compare_11_11( const FT& A1, const FT& B1,
	       const FT& A2, const FT& B2 )
{
  // Compares roots of (A1 X + B1) and (A2 X + B2).
  CGAL_precondition( A1 > 0 & A2 > 0 );
  return CGAL_NTS compare(B2*A1, B1*A2);
}

/*2 1*/template <class FT>
/*1 1*/Comparison_result
compare_21_11(const FT& A2, const FT& B2, const FT& C2,
              const FT& A1, const FT& B1 )
{
  // Compares roots of (A1 X + B1) and the smaller of (A2 X^2 + B2 X + C2).
  CGAL_precondition(A2 > 0);

  // First, we compare the root of P1 to the root of the derivative of P2.

  int cmp = compare_11_11<FT>(A1, B1, A2*2, B2);

  if (cmp > 0)
    return LARGER;

  // If it doesn't work, we evaluate the sign of P2 at the root of P1.

  FT p2 = B1 * (A1*B2 - A2*B1) - C2 * CGAL_NTS square(A1);

  return CGAL_NTS sign(p2);
}

/*2 2*/template <class FT>
/*2 1*/Comparison_result
compare_22_21( const FT& A1p, const FT& B1p, const FT& C1p,
	       const FT& A2p, const FT& B2p, const FT& C2p )
{
    // Compares the larger root of (A1 X^2 + B1 X + C1)
    //      to the smaller root of (A2 X^2 + B2 X + C2)
    // It boils down to the code from the DFMT paper
    // by multiplying A* and C* by 2, and B* by -1.

    CGAL_precondition(A1p > 0 & A2p > 0);

    FT A1 = 2 * A1p;
    FT C1 = 2 * C1p;
    FT B1 = -B1p;

    FT A2 = 2 * A2p;
    FT C2 = 2 * C2p;
    FT B2 = -B2p;

    // Now compares the larger root of (A1 X^2 -2B1 X + C1)
    //          to the smaller root of (A2 X^2 -2B2 X + C2)
    FT J = calcJ(A1,B1,A2,B2);

    if ( J < 0 ) return LARGER;   // r1 > l2

    FT K = calcK(A1,B1,C1,A2,B2,C2);

    if ( K < 0 ) return LARGER;   // r1 > l2

    FT Jp = calcJp(B1,C1,B2,C2);

    if ( Jp < 0 ) return SMALLER;  // r1 < l2

    FT P4 = calcP4(J,Jp,A1,C1,A2,C2);

    return - CGAL_NTS sign(P4);
    // if ( P4< FT(0) ) return LARGER;   // r1 > l2
    // if ( P4> FT(0) ) return SMALLER;  // r1 < l2
    // return EQUAL;
}

/*2 2*/template <class FT> inline
/*1 2*/Comparison_result
compare_22_12( const FT& A1, const FT& B1, const FT& C1,
	       const FT& A2, const FT& B2, const FT& C2 )
{
    // _22_12 boils down to _22_21 by :
    // - swapping the two polynomials
    // - changing the sign of the result
    return - compare_22_21(A2, B2, C2, A1, B1, C1);
}

/*2 2*/template <class FT>
/*1 1*/Comparison_result
compare_22_11( const FT& A1p, const FT& B1p, const FT& C1p,
	       const FT& A2p, const FT& B2p, const FT& C2p )
{
  // Compares the smaller root of (A1 X^2 + B1 X + C1)
  //       to the smaller root of (A2 X^2 + B2 X + C2)
  // It boils down to the code from the DFMT paper
  // by multiplying A* and C* by 2, and B* by -1.

  CGAL_precondition(A1p > 0 & A2p > 0);

  FT A1 = 2 * A1p;
  FT C1 = 2 * C1p;
  FT B1 = -B1p;

  FT A2 = 2 * A2p;
  FT C2 = 2 * C2p;
  FT B2 = -B2p;

  // Compares the smaller root of (A1 X^2 -2B1 X + C1)
  //       to the smaller root of (A2 X^2 -2B2 X + C2)
  FT J = calcJ(A1,B1,A2,B2);
  FT K = calcK(A1,B1,C1,A2,B2,C2);

  if (J > 0)
  {
    if (K > 0) return SMALLER;  // l1 < l2

    FT I1= calcI(A1,B1,C1);
    FT I2= calcI(A2,B2,C2);
    FT D = calcD(A1,I1,A2,I2);

    if (D > 0) return SMALLER;  // l1 < l2

    FT Jp = calcJp(B1,C1,B2,C2);

    if (Jp < 0) return LARGER;   // l1 > l2

    FT P4 = calcP4(I1,I2,K);

    return CGAL_NTS sign(P4);
  }

  // J <= 0
  if (K > 0) return LARGER;   // l1 > l2

  FT I1= calcI(A1,B1,C1);
  FT I2= calcI(A2,B2,C2);
  FT D = calcD(A1,I1,A2,I2);

  if (D < 0) return LARGER;   // l1 > l2

  FT Jp = calcJp(B1,C1,B2,C2);

  if (Jp > 0) return SMALLER;  // l1 < l2

  FT P4 = calcP4(I1,I2,K);

  return - CGAL_NTS sign(P4);
}

/*2 2*/template <class FT> inline
/*2 2*/Comparison_result
compare_22_22( const FT& A1, const FT& B1, const FT& C1,
	       const FT& A2, const FT& B2, const FT& C2 )
{
  // _22_22 boils down to _22_11 by :
  // - changing the sign of the two roots (X <-> -X in the polynomial)
  // - swapping the two polynomials
  return compare_22_11<FT>(A2, -B2, C2, A1, -B1, C1);
}

template <class FT>
inline FT calcI(const FT& A, const FT& B, const FT& C)
{ return CGAL_NTS square(B)-A*C; }

template <class FT>
inline FT calcJ(const FT& A1, const FT& B1, const FT& A2, const FT& B2)
{ return A1*B2-A2*B1; }

template <class FT>
inline FT calcK(const FT& A1, const FT& B1, const FT& C1,
		const FT& A2, const FT& B2, const FT& C2)
{ return C1*A2+A1*C2-2*B1*B2; }

template <class FT>
inline FT calcJp(const FT& B1, const FT& C1, const FT& B2, const FT& C2)
{ return B1*C2-C1*B2; }

template <class FT>
inline FT calcP4(const FT& J,  const FT& Jp,
		 const FT& A1, const FT& C1,
		 const FT& A2, const FT& C2)
{ return CGAL_NTS square(A1*C2-C1*A2)-4*J*Jp;}

template <class FT>
inline FT calcP4(const FT& I1, const FT& I2, const FT& K)
{ return CGAL_NTS square(K)-4*I1*I2;}

template <class FT>
inline FT calcD(const FT& A1, const FT& I1, const FT& A2, const FT& I2)
{ return I1*CGAL_NTS square(A2) - I2*CGAL_NTS square(A1);}

} // namespace internal
} // namespace CGAL

#endif // CGAL_NUMBER_TypeS_ROOT_OF_COMPARISON_FUNCTIONS_22_H
