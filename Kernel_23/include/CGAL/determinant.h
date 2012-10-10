// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// 
//
// Author(s)     : Sylvain Pion
//                 Stefan Schirra

#ifndef CGAL_DETERMINANT_H
#define CGAL_DETERMINANT_H

#include <CGAL/kernel_config.h>

namespace CGAL {

template <class RT>
inline
RT
determinant(
 const RT& a00,  const RT& a01,
 const RT& a10,  const RT& a11)
{
// First compute the det2x2
  const RT m01 = a00*a11 - a10*a01;
  return m01;
}

template <class RT>
CGAL_KERNEL_MEDIUM_INLINE
RT
determinant(
 const RT& a00,  const RT& a01,  const RT& a02,
 const RT& a10,  const RT& a11,  const RT& a12,
 const RT& a20,  const RT& a21,  const RT& a22)
{
// First compute the det2x2
  const RT m01 = a00*a11 - a10*a01;
  const RT m02 = a00*a21 - a20*a01;
  const RT m12 = a10*a21 - a20*a11;
// Now compute the minors of rank 3
  const RT m012 = m01*a22 - m02*a12 + m12*a02;
  return m012;
}

template <class RT>
CGAL_KERNEL_LARGE_INLINE
RT
determinant(
 const RT& a00,  const RT& a01,  const RT& a02,  const RT& a03,
 const RT& a10,  const RT& a11,  const RT& a12,  const RT& a13,
 const RT& a20,  const RT& a21,  const RT& a22,  const RT& a23,
 const RT& a30,  const RT& a31,  const RT& a32,  const RT& a33)
{
// First compute the det2x2
  const RT m01 = a10*a01 - a00*a11;
  const RT m02 = a20*a01 - a00*a21;
  const RT m03 = a30*a01 - a00*a31;
  const RT m12 = a20*a11 - a10*a21;
  const RT m13 = a30*a11 - a10*a31;
  const RT m23 = a30*a21 - a20*a31;
// Now compute the minors of rank 3
  const RT m012 = m12*a02 - m02*a12 + m01*a22;
  const RT m013 = m13*a02 - m03*a12 + m01*a32;
  const RT m023 = m23*a02 - m03*a22 + m02*a32;
  const RT m123 = m23*a12 - m13*a22 + m12*a32;
// Now compute the minors of rank 4
  const RT m0123 = m123*a03 - m023*a13 + m013*a23 - m012*a33;
  return m0123;
}

template <class RT>
CGAL_KERNEL_LARGE_INLINE
RT
determinant(
 const RT& a00,  const RT& a01,  const RT& a02,  const RT& a03,  const RT& a04,
 const RT& a10,  const RT& a11,  const RT& a12,  const RT& a13,  const RT& a14,
 const RT& a20,  const RT& a21,  const RT& a22,  const RT& a23,  const RT& a24,
 const RT& a30,  const RT& a31,  const RT& a32,  const RT& a33,  const RT& a34,
 const RT& a40,  const RT& a41,  const RT& a42,  const RT& a43,  const RT& a44)
{
// First compute the det2x2
  const RT m01 = a10*a01 - a00*a11;
  const RT m02 = a20*a01 - a00*a21;
  const RT m03 = a30*a01 - a00*a31;
  const RT m04 = a40*a01 - a00*a41;
  const RT m12 = a20*a11 - a10*a21;
  const RT m13 = a30*a11 - a10*a31;
  const RT m14 = a40*a11 - a10*a41;
  const RT m23 = a30*a21 - a20*a31;
  const RT m24 = a40*a21 - a20*a41;
  const RT m34 = a40*a31 - a30*a41;
// Now compute the minors of rank 3
  const RT m012 = m12*a02 - m02*a12 + m01*a22;
  const RT m013 = m13*a02 - m03*a12 + m01*a32;
  const RT m014 = m14*a02 - m04*a12 + m01*a42;
  const RT m023 = m23*a02 - m03*a22 + m02*a32;
  const RT m024 = m24*a02 - m04*a22 + m02*a42;
  const RT m034 = m34*a02 - m04*a32 + m03*a42;
  const RT m123 = m23*a12 - m13*a22 + m12*a32;
  const RT m124 = m24*a12 - m14*a22 + m12*a42;
  const RT m134 = m34*a12 - m14*a32 + m13*a42;
  const RT m234 = m34*a22 - m24*a32 + m23*a42;
// Now compute the minors of rank 4
  const RT m0123 = m123*a03 - m023*a13 + m013*a23 - m012*a33;
  const RT m0124 = m124*a03 - m024*a13 + m014*a23 - m012*a43;
  const RT m0134 = m134*a03 - m034*a13 + m014*a33 - m013*a43;
  const RT m0234 = m234*a03 - m034*a23 + m024*a33 - m023*a43;
  const RT m1234 = m234*a13 - m134*a23 + m124*a33 - m123*a43;
// Now compute the minors of rank 5
  const RT m01234 = m1234*a04 - m0234*a14 + m0134*a24 - m0124*a34 + m0123*a44;
  return m01234;
}

template <class RT>
RT
determinant(
 const RT& a00, const RT& a01, const RT& a02, const RT& a03, const RT& a04,
 const RT& a05,
 const RT& a10, const RT& a11, const RT& a12, const RT& a13, const RT& a14,
 const RT& a15,
 const RT& a20, const RT& a21, const RT& a22, const RT& a23, const RT& a24,
 const RT& a25,
 const RT& a30, const RT& a31, const RT& a32, const RT& a33, const RT& a34,
 const RT& a35,
 const RT& a40, const RT& a41, const RT& a42, const RT& a43, const RT& a44,
 const RT& a45,
 const RT& a50, const RT& a51, const RT& a52, const RT& a53, const RT& a54,
 const RT& a55)
{
// First compute the det2x2
  const RT m01 = a00*a11 - a10*a01;
  const RT m02 = a00*a21 - a20*a01;
  const RT m03 = a00*a31 - a30*a01;
  const RT m04 = a00*a41 - a40*a01;
  const RT m05 = a00*a51 - a50*a01;
  const RT m12 = a10*a21 - a20*a11;
  const RT m13 = a10*a31 - a30*a11;
  const RT m14 = a10*a41 - a40*a11;
  const RT m15 = a10*a51 - a50*a11;
  const RT m23 = a20*a31 - a30*a21;
  const RT m24 = a20*a41 - a40*a21;
  const RT m25 = a20*a51 - a50*a21;
  const RT m34 = a30*a41 - a40*a31;
  const RT m35 = a30*a51 - a50*a31;
  const RT m45 = a40*a51 - a50*a41;
// Now compute the minors of rank 3
  const RT m012 = m01*a22 - m02*a12 + m12*a02;
  const RT m013 = m01*a32 - m03*a12 + m13*a02;
  const RT m014 = m01*a42 - m04*a12 + m14*a02;
  const RT m015 = m01*a52 - m05*a12 + m15*a02;
  const RT m023 = m02*a32 - m03*a22 + m23*a02;
  const RT m024 = m02*a42 - m04*a22 + m24*a02;
  const RT m025 = m02*a52 - m05*a22 + m25*a02;
  const RT m034 = m03*a42 - m04*a32 + m34*a02;
  const RT m035 = m03*a52 - m05*a32 + m35*a02;
  const RT m045 = m04*a52 - m05*a42 + m45*a02;
  const RT m123 = m12*a32 - m13*a22 + m23*a12;
  const RT m124 = m12*a42 - m14*a22 + m24*a12;
  const RT m125 = m12*a52 - m15*a22 + m25*a12;
  const RT m134 = m13*a42 - m14*a32 + m34*a12;
  const RT m135 = m13*a52 - m15*a32 + m35*a12;
  const RT m145 = m14*a52 - m15*a42 + m45*a12;
  const RT m234 = m23*a42 - m24*a32 + m34*a22;
  const RT m235 = m23*a52 - m25*a32 + m35*a22;
  const RT m245 = m24*a52 - m25*a42 + m45*a22;
  const RT m345 = m34*a52 - m35*a42 + m45*a32;
// Now compute the minors of rank 4
  const RT m0123 = m012*a33 - m013*a23 + m023*a13 - m123*a03;
  const RT m0124 = m012*a43 - m014*a23 + m024*a13 - m124*a03;
  const RT m0125 = m012*a53 - m015*a23 + m025*a13 - m125*a03;
  const RT m0134 = m013*a43 - m014*a33 + m034*a13 - m134*a03;
  const RT m0135 = m013*a53 - m015*a33 + m035*a13 - m135*a03;
  const RT m0145 = m014*a53 - m015*a43 + m045*a13 - m145*a03;
  const RT m0234 = m023*a43 - m024*a33 + m034*a23 - m234*a03;
  const RT m0235 = m023*a53 - m025*a33 + m035*a23 - m235*a03;
  const RT m0245 = m024*a53 - m025*a43 + m045*a23 - m245*a03;
  const RT m0345 = m034*a53 - m035*a43 + m045*a33 - m345*a03;
  const RT m1234 = m123*a43 - m124*a33 + m134*a23 - m234*a13;
  const RT m1235 = m123*a53 - m125*a33 + m135*a23 - m235*a13;
  const RT m1245 = m124*a53 - m125*a43 + m145*a23 - m245*a13;
  const RT m1345 = m134*a53 - m135*a43 + m145*a33 - m345*a13;
  const RT m2345 = m234*a53 - m235*a43 + m245*a33 - m345*a23;
// Now compute the minors of rank 5
  const RT m01234 = m0123*a44 - m0124*a34 + m0134*a24 - m0234*a14 + m1234*a04;
  const RT m01235 = m0123*a54 - m0125*a34 + m0135*a24 - m0235*a14 + m1235*a04;
  const RT m01245 = m0124*a54 - m0125*a44 + m0145*a24 - m0245*a14 + m1245*a04;
  const RT m01345 = m0134*a54 - m0135*a44 + m0145*a34 - m0345*a14 + m1345*a04;
  const RT m02345 = m0234*a54 - m0235*a44 + m0245*a34 - m0345*a24 + m2345*a04;
  const RT m12345 = m1234*a54 - m1235*a44 + m1245*a34 - m1345*a24 + m2345*a14;
// Now compute the minors of rank 6
  const RT m012345 = m01234*a55 - m01235*a45 + m01245*a35 - m01345*a25
                   + m02345*a15 - m12345*a05;
  return m012345;
}

} //namespace CGAL

#endif // CGAL_DETERMINANT_H
