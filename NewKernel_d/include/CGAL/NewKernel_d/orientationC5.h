// Copyright (c) 2025
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Andreas Fabri

// Files needed as input for the static filter generator

#if 0

double
determinant(
 double a00,  double a01,  double a02,  double a03,  double a04,
 double a10,  double a11,  double a12,  double a13,  double a14,
 double a20,  double a21,  double a22,  double a23,  double a24,
 double a30,  double a31,  double a32,  double a33,  double a34,
 double a40,  double a41,  double a42,  double a43,  double a44)
{
  double m01 = a10*a01 - a00*a11;
  double m02 = a20*a01 - a00*a21;
  double m03 = a30*a01 - a00*a31;
  double m04 = a40*a01 - a00*a41;
  double m12 = a20*a11 - a10*a21;
  double m13 = a30*a11 - a10*a31;
  double m14 = a40*a11 - a10*a41;
  double m23 = a30*a21 - a20*a31;
  double m24 = a40*a21 - a20*a41;
  double m34 = a40*a31 - a30*a41;
  double m012 = m12*a02 - m02*a12 + m01*a22;
  double m013 = m13*a02 - m03*a12 + m01*a32;
  double m014 = m14*a02 - m04*a12 + m01*a42;
  double m023 = m23*a02 - m03*a22 + m02*a32;
  double m024 = m24*a02 - m04*a22 + m02*a42;
  double m034 = m34*a02 - m04*a32 + m03*a42;
  double m123 = m23*a12 - m13*a22 + m12*a32;
  double m124 = m24*a12 - m14*a22 + m12*a42;
  double m134 = m34*a12 - m14*a32 + m13*a42;
  double m234 = m34*a22 - m24*a32 + m23*a42;
  double m0123 = m123*a03 - m023*a13 + m013*a23 - m012*a33;
  double m0124 = m124*a03 - m024*a13 + m014*a23 - m012*a43;
  double m0134 = m134*a03 - m034*a13 + m014*a33 - m013*a43;
  double m0234 = m234*a03 - m034*a23 + m024*a33 - m023*a43;
  double m1234 = m234*a13 - m134*a23 + m124*a33 - m123*a43;
  double m01234 = m1234*a04 - m0234*a14 + m0134*a24 - m0124*a34 + m0123*a44;
  return m01234;
}

int orientationC5(double p0, double p1, double p2, double p3, double p4,
                  double q0, double q1, double q2, double q3, double q4,
                  double r0, double r1, double r2, double r3, double r4,
                  double s0, double s1, double s2, double s3, double s4,
                  double t0, double t1, double t2, double t3, double t4,
                  double u0, double u1, double u2, double u3, double u4)
group p0 q0 r0 t0 q0 u0;
group p1 q1 r1 t1 q1 u1;
group p2 q2 r2 t2 q2 u2;
group p3 q3 r3 t3 q3 u3;
group p4 q4 r4 t4 q4 u4;
{
    double m01;
    m01 = (q0 - p0);
    double m02;
    m02 = (r0 - p0);
    double m03;
    m03 = (s0 - p0);
    double m04;
    m04 = (t0 - p0);
    double m05;
    m05 = (u0 - p0);

    double m11;
    m11 = (q1 - p1);
    double m12;
    m12 = (r1 - p1);
    double m13;
    m13 = (s1 - p1);
    double m14;
    m14 = (t1 - p1);
    double m15;
    m15 = (u1 - p1);

    double m21;
    m21 = (q2 - p2);
    double m22;
    m22 = (r2 - p2);
    double m23;
    m23 = (s2 - p2);
    double m24;
    m24 = (t2 - p2);
    double m25;
    m25 = (u2 - p2);

    double m31;
    m31 = (q3 - p3);
    double m32;
    m32 = (r3 - p3);
    double m33;
    m33 = (s3 - p3);
    double m34;
    m34 = (t3 - p3);
    double m35;
    m35 = (u3 - p3);

    double m41;
    m41 = (q4 - p4);
    double m42;
    m42 = (r4 - p4);
    double m43;
    m43 = (s4 - p4);
    double m44;
    m44 = (t4 - p4);
    double m45;
    m45 = (u4 - p4);

    double det = determinant(m01, m02, m03, m04, m05,
                             m11, m12, m13, m14, m15,
                             m21, m22, m23, m24, m25,
                             m31, m32, m33, m34, m35,
                             m41, m42, m43, m44, m45);
    return sign(det);
}

//===========generated ==================

inline int orientationC5( double p0, double p1, double p2, double p3, double p4, double q0, double q1, double q2, double q3, double q4, double r0, double r1, double r2, double r3, double r4, double s0, double s1, double s2, double s3, double s4, double t0, double t1, double t2, double t3, double t4, double u0, double u1, double u2, double u3, double u4) {
    double m01;
    m01 = (q0 - p0);
    double m02;
    m02 = (r0 - p0);
    double m03;
    m03 = (s0 - p0);
    double m04;
    m04 = (t0 - p0);
    double m05;
    m05 = (u0 - p0);
    double m11;
    m11 = (q1 - p1);
    double m12;
    m12 = (r1 - p1);
    double m13;
    m13 = (s1 - p1);
    double m14;
    m14 = (t1 - p1);
    double m15;
    m15 = (u1 - p1);
    double m21;
    m21 = (q2 - p2);
    double m22;
    m22 = (r2 - p2);
    double m23;
    m23 = (s2 - p2);
    double m24;
    m24 = (t2 - p2);
    double m25;
    m25 = (u2 - p2);
    double m31;
    m31 = (q3 - p3);
    double m32;
    m32 = (r3 - p3);
    double m33;
    m33 = (s3 - p3);
    double m34;
    m34 = (t3 - p3);
    double m35;
    m35 = (u3 - p3);
    double m41;
    m41 = (q4 - p4);
    double m42;
    m42 = (r4 - p4);
    double m43;
    m43 = (s4 - p4);
    double m44;
    m44 = (t4 - p4);
    double m45;
    m45 = (u4 - p4);
    double det;
    det = determinant( m01, m02, m03, m04, m05, m11, m12, m13, m14, m15, m21, m22, m23, m24, m25, m31, m32, m33, m34, m35, m41, m42, m43, m44, m45 );
    int int_tmp_result;
    double eps;
    double max1 = CGAL::abs(m01);
    double am = CGAL::abs(m02);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m03);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m11);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m12);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m13);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m21);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m22);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m23);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m33);
    if( (max1 < am) ) { max1 = am; }


    double max2 = CGAL::abs(m01);
    am = CGAL::abs(m02);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m11);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m12);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m21);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m22);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m23);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m31);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m32);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m33);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m43);
    if( (max2 < am) ) { max2 = am; }


    double max3 = CGAL::abs(m04);
    am = CGAL::abs(m14);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m24);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m34);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m44);
    if( (max3 < am) ) { max3 = am; }


    double max4 = CGAL::abs(m05);
    am = CGAL::abs(m15);
    if( (max4 < am) ) { max4 = am; }
    am = CGAL::abs(m25);
    if( (max4 < am) ) { max4 = am; }
    am = CGAL::abs(m35);
    if( (max4 < am) ) { max4 = am; }
    am = CGAL::abs(m45);
    if( (max4 < am) ) { max4 = am; }


    double max5 = CGAL::abs(m11);
    am = CGAL::abs(m12);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m21);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m22);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m31);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m32);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m41);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m42);
    if( (max5 < am) ) { max5 = am; }

    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max5;
    upper_bound_1 = max5;
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (lower_bound_1 < 9.99657131447050328602e-60) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 3.21387608851797912384e+60) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (2.22889232457534740153e-13 * ((((max5 * max2) * max1) * max3) * max4));
        if( (det > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (det < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}

#endif
