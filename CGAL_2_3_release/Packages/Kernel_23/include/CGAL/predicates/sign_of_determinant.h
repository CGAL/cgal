// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : predicates/sign_of_determinant.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_PREDICATES_SIGN_OF_DETERMINANT_H
#define CGAL_PREDICATES_SIGN_OF_DETERMINANT_H

#include <CGAL/determinant.h>

CGAL_BEGIN_NAMESPACE

template <class FT>
inline
Sign
sign_of_determinant2x2( const FT& a00,  const FT& a01,
                        const FT& a10,  const FT& a11)
{
  return
    static_cast<Sign>(static_cast<int>(CGAL_NTS compare( a00*a11, a10*a01)));
}

template <class FT>
inline
Sign
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
Sign
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
Sign
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
Sign
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
