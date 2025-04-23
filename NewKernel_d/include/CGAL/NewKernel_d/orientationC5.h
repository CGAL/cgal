


typedef double RT;

RT
determinant(
 RT a00,  RT a01,  RT a02,  RT a03,  RT a04,
 RT a10,  RT a11,  RT a12,  RT a13,  RT a14,
 RT a20,  RT a21,  RT a22,  RT a23,  RT a24,
 RT a30,  RT a31,  RT a32,  RT a33,  RT a34,
 RT a40,  RT a41,  RT a42,  RT a43,  RT a44)
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



int
orientationC5(RT p0, RT p1, RT p2, RT p3, RT p4,
              RT q0, RT q1, RT q2, RT q3, RT q4,
              RT r0, RT r1, RT r2, RT r3, RT r4,
              RT s0, RT s1, RT s2, RT s3, RT s4,
              RT t0, RT t1, RT t2, RT t3, RT t4,
              RT u0, RT u1, RT u2, RT u3, RT u4)
{
  RT pq0 = q0 - p0;
  RT pq1 = q1 - p1;
  RT pq2 = q2 - p2;
  RT pq3 = q3 - p3;
  RT pq4 = q4 - p4;
  RT pr0 = r0 - p0;
  RT pr1 = r1 - p1;
  RT pr2 = r2 - p2;
  RT pr3 = r3 - p3;
  RT pr4 = r4 - p4;
  RT ps0 = s0 - p0;
  RT ps1 = s1 - p1;
  RT ps2 = s2 - p2;
  RT ps3 = s3 - p3;
  RT ps4 = s4 - p4;
  RT pt0 = t0 - p0;
  RT pt1 = t1 - p1;
  RT pt2 = t2 - p2;
  RT pt3 = t3 - p3;
  RT pt4 = t4 - p4;
  RT pu0 = t0 - p0;
  RT pu1 = u1 - p1;
  RT pu2 = u2 - p2;
  RT pu3 = u3 - p3;
  RT pu4 = u4 - p4;
  RT det = determinant(pq0, pq1, pq2, pq3, pq4,
                       pr0, pr1, pr2, pr3, pr4,
                       ps0, ps1, ps2, ps3, ps4,
                       pt0, pt1, pt2, pt3, pt4,
                       pu0, pu1, pu2, pu3, pu4);

  if (det > 0) return 1;
  if (det < 0) return -1;
}
