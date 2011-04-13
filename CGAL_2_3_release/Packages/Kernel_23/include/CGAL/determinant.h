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
// file          : determinant.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_DETERMINANT_H
#define CGAL_DETERMINANT_H

CGAL_BEGIN_NAMESPACE

template <class FT>
inline
FT
det2x2_by_formula(
 const FT& a00,  const FT& a01,
 const FT& a10,  const FT& a11)
{
// First compute the det2x2
  const FT m01 = a00*a11 - a10*a01;
  return m01;
}

template <class FT>
CGAL_KERNEL_MEDIUM_INLINE
FT
det3x3_by_formula(
 const FT& a00,  const FT& a01,  const FT& a02,
 const FT& a10,  const FT& a11,  const FT& a12,
 const FT& a20,  const FT& a21,  const FT& a22)
{
// First compute the det2x2
  const FT m01 = a00*a11 - a10*a01;
  const FT m02 = a00*a21 - a20*a01;
  const FT m12 = a10*a21 - a20*a11;
// Now compute the minors of rank 3
  const FT m012 = m01*a22 - m02*a12 + m12*a02;
  return m012;
}

template <class FT>
CGAL_KERNEL_LARGE_INLINE
FT
det4x4_by_formula(
 const FT& a00,  const FT& a01,  const FT& a02,  const FT& a03,
 const FT& a10,  const FT& a11,  const FT& a12,  const FT& a13,
 const FT& a20,  const FT& a21,  const FT& a22,  const FT& a23,
 const FT& a30,  const FT& a31,  const FT& a32,  const FT& a33)
{
// First compute the det2x2
  const FT m01 = a10*a01 - a00*a11;
  const FT m02 = a20*a01 - a00*a21;
  const FT m03 = a30*a01 - a00*a31;
  const FT m12 = a20*a11 - a10*a21;
  const FT m13 = a30*a11 - a10*a31;
  const FT m23 = a30*a21 - a20*a31;
// Now compute the minors of rank 3
  const FT m012 = m12*a02 - m02*a12 + m01*a22;
  const FT m013 = m13*a02 - m03*a12 + m01*a32;
  const FT m023 = m23*a02 - m03*a22 + m02*a32;
  const FT m123 = m23*a12 - m13*a22 + m12*a32;
// Now compute the minors of rank 4
  const FT m0123 = m123*a03 - m023*a13 + m013*a23 - m012*a33;
  return m0123;
}

template <class FT>
CGAL_KERNEL_LARGE_INLINE
FT
det5x5_by_formula(
 const FT& a00,  const FT& a01,  const FT& a02,  const FT& a03,  const FT& a04,
 const FT& a10,  const FT& a11,  const FT& a12,  const FT& a13,  const FT& a14,
 const FT& a20,  const FT& a21,  const FT& a22,  const FT& a23,  const FT& a24,
 const FT& a30,  const FT& a31,  const FT& a32,  const FT& a33,  const FT& a34,
 const FT& a40,  const FT& a41,  const FT& a42,  const FT& a43,  const FT& a44)
{
// First compute the det2x2
  const FT m01 = a10*a01 - a00*a11;
  const FT m02 = a20*a01 - a00*a21;
  const FT m03 = a30*a01 - a00*a31;
  const FT m04 = a40*a01 - a00*a41;
  const FT m12 = a20*a11 - a10*a21;
  const FT m13 = a30*a11 - a10*a31;
  const FT m14 = a40*a11 - a10*a41;
  const FT m23 = a30*a21 - a20*a31;
  const FT m24 = a40*a21 - a20*a41;
  const FT m34 = a40*a31 - a30*a41;
// Now compute the minors of rank 3
  const FT m012 = m12*a02 - m02*a12 + m01*a22;
  const FT m013 = m13*a02 - m03*a12 + m01*a32;
  const FT m014 = m14*a02 - m04*a12 + m01*a42;
  const FT m023 = m23*a02 - m03*a22 + m02*a32;
  const FT m024 = m24*a02 - m04*a22 + m02*a42;
  const FT m034 = m34*a02 - m04*a32 + m03*a42;
  const FT m123 = m23*a12 - m13*a22 + m12*a32;
  const FT m124 = m24*a12 - m14*a22 + m12*a42;
  const FT m134 = m34*a12 - m14*a32 + m13*a42;
  const FT m234 = m34*a22 - m24*a32 + m23*a42;
// Now compute the minors of rank 4
  const FT m0123 = m123*a03 - m023*a13 + m013*a23 - m012*a33;
  const FT m0124 = m124*a03 - m024*a13 + m014*a23 - m012*a43;
  const FT m0134 = m134*a03 - m034*a13 + m014*a33 - m013*a43;
  const FT m0234 = m234*a03 - m034*a23 + m024*a33 - m023*a43;
  const FT m1234 = m234*a13 - m134*a23 + m124*a33 - m123*a43;
// Now compute the minors of rank 5
  const FT m01234 = m1234*a04 - m0234*a14 + m0134*a24 - m0124*a34 + m0123*a44;
  return m01234;
}

template <class FT>
FT
det6x6_by_formula(
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
// First compute the det2x2
  const FT m01 = a00*a11 - a10*a01;
  const FT m02 = a00*a21 - a20*a01;
  const FT m03 = a00*a31 - a30*a01;
  const FT m04 = a00*a41 - a40*a01;
  const FT m05 = a00*a51 - a50*a01;
  const FT m12 = a10*a21 - a20*a11;
  const FT m13 = a10*a31 - a30*a11;
  const FT m14 = a10*a41 - a40*a11;
  const FT m15 = a10*a51 - a50*a11;
  const FT m23 = a20*a31 - a30*a21;
  const FT m24 = a20*a41 - a40*a21;
  const FT m25 = a20*a51 - a50*a21;
  const FT m34 = a30*a41 - a40*a31;
  const FT m35 = a30*a51 - a50*a31;
  const FT m45 = a40*a51 - a50*a41;
// Now compute the minors of rank 3
  const FT m012 = m01*a22 - m02*a12 + m12*a02;
  const FT m013 = m01*a32 - m03*a12 + m13*a02;
  const FT m014 = m01*a42 - m04*a12 + m14*a02;
  const FT m015 = m01*a52 - m05*a12 + m15*a02;
  const FT m023 = m02*a32 - m03*a22 + m23*a02;
  const FT m024 = m02*a42 - m04*a22 + m24*a02;
  const FT m025 = m02*a52 - m05*a22 + m25*a02;
  const FT m034 = m03*a42 - m04*a32 + m34*a02;
  const FT m035 = m03*a52 - m05*a32 + m35*a02;
  const FT m045 = m04*a52 - m05*a42 + m45*a02;
  const FT m123 = m12*a32 - m13*a22 + m23*a12;
  const FT m124 = m12*a42 - m14*a22 + m24*a12;
  const FT m125 = m12*a52 - m15*a22 + m25*a12;
  const FT m134 = m13*a42 - m14*a32 + m34*a12;
  const FT m135 = m13*a52 - m15*a32 + m35*a12;
  const FT m145 = m14*a52 - m15*a42 + m45*a12;
  const FT m234 = m23*a42 - m24*a32 + m34*a22;
  const FT m235 = m23*a52 - m25*a32 + m35*a22;
  const FT m245 = m24*a52 - m25*a42 + m45*a22;
  const FT m345 = m34*a52 - m35*a42 + m45*a32;
// Now compute the minors of rank 4
  const FT m0123 = m012*a33 - m013*a23 + m023*a13 - m123*a03;
  const FT m0124 = m012*a43 - m014*a23 + m024*a13 - m124*a03;
  const FT m0125 = m012*a53 - m015*a23 + m025*a13 - m125*a03;
  const FT m0134 = m013*a43 - m014*a33 + m034*a13 - m134*a03;
  const FT m0135 = m013*a53 - m015*a33 + m035*a13 - m135*a03;
  const FT m0145 = m014*a53 - m015*a43 + m045*a13 - m145*a03;
  const FT m0234 = m023*a43 - m024*a33 + m034*a23 - m234*a03;
  const FT m0235 = m023*a53 - m025*a33 + m035*a23 - m235*a03;
  const FT m0245 = m024*a53 - m025*a43 + m045*a23 - m245*a03;
  const FT m0345 = m034*a53 - m035*a43 + m045*a33 - m345*a03;
  const FT m1234 = m123*a43 - m124*a33 + m134*a23 - m234*a13;
  const FT m1235 = m123*a53 - m125*a33 + m135*a23 - m235*a13;
  const FT m1245 = m124*a53 - m125*a43 + m145*a23 - m245*a13;
  const FT m1345 = m134*a53 - m135*a43 + m145*a33 - m345*a13;
  const FT m2345 = m234*a53 - m235*a43 + m245*a33 - m345*a23;
// Now compute the minors of rank 5
  const FT m01234 = m0123*a44 - m0124*a34 + m0134*a24 - m0234*a14 + m1234*a04;
  const FT m01235 = m0123*a54 - m0125*a34 + m0135*a24 - m0235*a14 + m1235*a04;
  const FT m01245 = m0124*a54 - m0125*a44 + m0145*a24 - m0245*a14 + m1245*a04;
  const FT m01345 = m0134*a54 - m0135*a44 + m0145*a34 - m0345*a14 + m1345*a04;
  const FT m02345 = m0234*a54 - m0235*a44 + m0245*a34 - m0345*a24 + m2345*a04;
  const FT m12345 = m1234*a54 - m1235*a44 + m1245*a34 - m1345*a24 + m2345*a14;
// Now compute the minors of rank 6
  const FT m012345 = m01234*a55 - m01235*a45 + m01245*a35 - m01345*a25
                   + m02345*a15 - m12345*a05;
  return m012345;
}

CGAL_END_NAMESPACE

#endif // CGAL_DETERMINANT_H
