// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion
//                 Stefan Schirra

#ifndef CGAL_PREDICATES_SIGN_OF_DETERMINANT_H
#define CGAL_PREDICATES_SIGN_OF_DETERMINANT_H

#include <CGAL/determinant.h>

CGAL_BEGIN_NAMESPACE

template <class FT>
inline
typename Sgn<FT>::result_type
sign_of_determinant2x2( const FT& a00,  const FT& a01,
                        const FT& a10,  const FT& a11)
{
  return enum_cast<Sign>(CGAL_NTS compare( a00*a11, a10*a01));
}

template <class FT>
inline
typename Sgn<FT>::result_type
sign_of_determinant3x3( const FT& a00,  const FT& a01,  const FT& a02,
                        const FT& a10,  const FT& a11,  const FT& a12,
                        const FT& a20,  const FT& a21,  const FT& a22)
{
  return CGAL_NTS sign(det3x3_by_formula(a00, a01, a02,
                                         a10, a11, a12,
                                         a20, a21, a22));
}

template <class FT>
inline
typename Sgn<FT>::result_type
sign_of_determinant4x4(
 const FT& a00,  const FT& a01,  const FT& a02,  const FT& a03,
 const FT& a10,  const FT& a11,  const FT& a12,  const FT& a13,
 const FT& a20,  const FT& a21,  const FT& a22,  const FT& a23,
 const FT& a30,  const FT& a31,  const FT& a32,  const FT& a33)
{
  return CGAL_NTS sign(det4x4_by_formula(a00, a01, a02, a03,
                                         a10, a11, a12, a13,
                                         a20, a21, a22, a23,
                                         a30, a31, a32, a33));
}

template <class FT>
CGAL_KERNEL_LARGE_INLINE
typename Sgn<FT>::result_type
sign_of_determinant5x5(
 const FT& a00,  const FT& a01,  const FT& a02,  const FT& a03,  const FT& a04,
 const FT& a10,  const FT& a11,  const FT& a12,  const FT& a13,  const FT& a14,
 const FT& a20,  const FT& a21,  const FT& a22,  const FT& a23,  const FT& a24,
 const FT& a30,  const FT& a31,  const FT& a32,  const FT& a33,  const FT& a34,
 const FT& a40,  const FT& a41,  const FT& a42,  const FT& a43,  const FT& a44)
{
  return CGAL_NTS sign(det5x5_by_formula(a00, a01, a02, a03, a04,
                                         a10, a11, a12, a13, a14,
                                         a20, a21, a22, a23, a24,
                                         a30, a31, a32, a33, a34,
                                         a40, a41, a42, a43, a44));
}

template <class FT>
CGAL_KERNEL_LARGE_INLINE
typename Sgn<FT>::result_type
sign_of_determinant6x6(
 const FT& a00, const FT& a01, const FT& a02, const FT& a03, const FT& a04,
 const FT& a05,
 const FT& a10, const FT& a11, const FT& a12, const FT& a13, const FT& a14,
 const FT& a15,
 const FT& a20, const FT& a21, const FT& a22, const FT& a23, const FT& a24,
 const FT& a25,
 const FT& a30, const FT& a31, const FT& a32, const FT& a33, const FT& a34,
 const FT& a35,
 const FT& a40, const FT& a41, const FT& a42, const FT& a43, const FT& a44,
 const FT& a45,
 const FT& a50, const FT& a51, const FT& a52, const FT& a53, const FT& a54,
 const FT& a55)
{
  return CGAL_NTS sign(det6x6_by_formula(a00, a01, a02, a03, a04, a05,
                                         a10, a11, a12, a13, a14, a15,
                                         a20, a21, a22, a23, a24, a25,
                                         a30, a31, a32, a33, a34, a35,
                                         a40, a41, a42, a43, a44, a45,
                                         a50, a51, a52, a53, a54, a55));
}

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#include <CGAL/Arithmetic_filter/predicates/sign_of_determinant.h>
#endif // CGAL_ARITHMETIC_FILTER_H

#endif // CGAL_PREDICATES_SIGN_OF_DETERMINANT_H
