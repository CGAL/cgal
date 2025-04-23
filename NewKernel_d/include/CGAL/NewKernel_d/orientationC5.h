// Copyright (c) 2025  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

typedef double RT;

RT
determinant_with_1_in_row_0(
 RT a01, RT a02, RT a03, RT a04, RT a05,
 RT a11, RT a12, RT a13, RT a14, RT a15,
 RT a21, RT a22, RT a23, RT a24, RT a25,
 RT a31, RT a32, RT a33, RT a34, RT a35,
 RT a41, RT a42, RT a43, RT a44, RT a45,
 RT a51, RT a52, RT a53, RT a54, RT a55)
{
// First compute the det2x2
  const RT m01 = a11 - a01;
  const RT m02 = a21 - a01;
  const RT m03 = a31 - a01;
  const RT m04 = a41 - a01;
  const RT m05 = a51 - a01;
  const RT m12 = a21 - a11;
  const RT m13 = a31 - a11;
  const RT m14 = a41 - a11;
  const RT m15 = a51 - a11;
  const RT m23 = a31 - a21;
  const RT m24 = a41 - a21;
  const RT m25 = a51 - a21;
  const RT m34 = a41 - a31;
  const RT m35 = a51 - a31;
  const RT m45 = a51 - a41;
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



inline int orientationC5( double p0, double p1, double p2, double p3, double p4,
                          double q0, double q1, double q2, double q3, double q4,
                          double r0, double r1, double r2, double r3, double r4,
                          double s0, double s1, double s2, double s3, double s4,
                          double t0, double t1, double t2, double t3, double t4,
                          double u0, double u1, double u2, double u3, double u4)
{
  RT det = determinant_with_1_in_row_0(p0, p1, p2, p3, p4,
                                       q0, q1, q2, q3, q4,
                                       r0, r1, r2, r3, r4,
                                       s0, s1, s2, s3, s4,
                                       t0, t1, t2, t3, t4,
                                       u0, u1, u2, u3, u4);
    if (det > 0) return 1;
    if (det < 0) return -1;
    return 0;
}

