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


inline double determinant_with_1_in_row_0( double a01, double a02, double a03, double a04, double a05, double a11, double a12, double a13, double a14, double a15, double a21, double a22, double a23, double a24, double a25, double a31, double a32, double a33, double a34, double a35, double a41, double a42, double a43, double a44, double a45, double a51, double a52, double a53, double a54, double a55) {
    double m01;
    m01 = (a11 - a01);
    double m02;
    m02 = (a21 - a01);
    double m03;
    m03 = (a31 - a01);
    double m04;
    m04 = (a41 - a01);
    double m05;
    m05 = (a51 - a01);
    double m12;
    m12 = (a21 - a11);
    double m13;
    m13 = (a31 - a11);
    double m14;
    m14 = (a41 - a11);
    double m15;
    m15 = (a51 - a11);
    double m23;
    m23 = (a31 - a21);
    double m24;
    m24 = (a41 - a21);
    double m25;
    m25 = (a51 - a21);
    double m34;
    m34 = (a41 - a31);
    double m35;
    m35 = (a51 - a31);
    double m45;
    m45 = (a51 - a41);
    double m012;
    m012 = (((m01 * a22) - (m02 * a12)) + (m12 * a02));
    double m013;
    m013 = (((m01 * a32) - (m03 * a12)) + (m13 * a02));
    double m014;
    m014 = (((m01 * a42) - (m04 * a12)) + (m14 * a02));
    double m015;
    m015 = (((m01 * a52) - (m05 * a12)) + (m15 * a02));
    double m023;
    m023 = (((m02 * a32) - (m03 * a22)) + (m23 * a02));
    double m024;
    m024 = (((m02 * a42) - (m04 * a22)) + (m24 * a02));
    double m025;
    m025 = (((m02 * a52) - (m05 * a22)) + (m25 * a02));
    double m034;
    m034 = (((m03 * a42) - (m04 * a32)) + (m34 * a02));
    double m035;
    m035 = (((m03 * a52) - (m05 * a32)) + (m35 * a02));
    double m045;
    m045 = (((m04 * a52) - (m05 * a42)) + (m45 * a02));
    double m123;
    m123 = (((m12 * a32) - (m13 * a22)) + (m23 * a12));
    double m124;
    m124 = (((m12 * a42) - (m14 * a22)) + (m24 * a12));
    double m125;
    m125 = (((m12 * a52) - (m15 * a22)) + (m25 * a12));
    double m134;
    m134 = (((m13 * a42) - (m14 * a32)) + (m34 * a12));
    double m135;
    m135 = (((m13 * a52) - (m15 * a32)) + (m35 * a12));
    double m145;
    m145 = (((m14 * a52) - (m15 * a42)) + (m45 * a12));
    double m234;
    m234 = (((m23 * a42) - (m24 * a32)) + (m34 * a22));
    double m235;
    m235 = (((m23 * a52) - (m25 * a32)) + (m35 * a22));
    double m245;
    m245 = (((m24 * a52) - (m25 * a42)) + (m45 * a22));
    double m345;
    m345 = (((m34 * a52) - (m35 * a42)) + (m45 * a32));
    double m0123;
    m0123 = ((((m012 * a33) - (m013 * a23)) + (m023 * a13)) - (m123 * a03));
    double m0124;
    m0124 = ((((m012 * a43) - (m014 * a23)) + (m024 * a13)) - (m124 * a03));
    double m0125;
    m0125 = ((((m012 * a53) - (m015 * a23)) + (m025 * a13)) - (m125 * a03));
    double m0134;
    m0134 = ((((m013 * a43) - (m014 * a33)) + (m034 * a13)) - (m134 * a03));
    double m0135;
    m0135 = ((((m013 * a53) - (m015 * a33)) + (m035 * a13)) - (m135 * a03));
    double m0145;
    m0145 = ((((m014 * a53) - (m015 * a43)) + (m045 * a13)) - (m145 * a03));
    double m0234;
    m0234 = ((((m023 * a43) - (m024 * a33)) + (m034 * a23)) - (m234 * a03));
    double m0235;
    m0235 = ((((m023 * a53) - (m025 * a33)) + (m035 * a23)) - (m235 * a03));
    double m0245;
    m0245 = ((((m024 * a53) - (m025 * a43)) + (m045 * a23)) - (m245 * a03));
    double m0345;
    m0345 = ((((m034 * a53) - (m035 * a43)) + (m045 * a33)) - (m345 * a03));
    double m1234;
    m1234 = ((((m123 * a43) - (m124 * a33)) + (m134 * a23)) - (m234 * a13));
    double m1235;
    m1235 = ((((m123 * a53) - (m125 * a33)) + (m135 * a23)) - (m235 * a13));
    double m1245;
    m1245 = ((((m124 * a53) - (m125 * a43)) + (m145 * a23)) - (m245 * a13));
    double m1345;
    m1345 = ((((m134 * a53) - (m135 * a43)) + (m145 * a33)) - (m345 * a13));
    double m2345;
    m2345 = ((((m234 * a53) - (m235 * a43)) + (m245 * a33)) - (m345 * a23));
    double m01234;
    m01234 = (((((m0123 * a44) - (m0124 * a34)) + (m0134 * a24)) - (m0234 * a14)) + (m1234 * a04));
    double m01235;
    m01235 = (((((m0123 * a54) - (m0125 * a34)) + (m0135 * a24)) - (m0235 * a14)) + (m1235 * a04));
    double m01245;
    m01245 = (((((m0124 * a54) - (m0125 * a44)) + (m0145 * a24)) - (m0245 * a14)) + (m1245 * a04));
    double m01345;
    m01345 = (((((m0134 * a54) - (m0135 * a44)) + (m0145 * a34)) - (m0345 * a14)) + (m1345 * a04));
    double m02345;
    m02345 = (((((m0234 * a54) - (m0235 * a44)) + (m0245 * a34)) - (m0345 * a24)) + (m2345 * a04));
    double m12345;
    m12345 = (((((m1234 * a54) - (m1235 * a44)) + (m1245 * a34)) - (m1345 * a24)) + (m2345 * a14));
    double m012345;
    m012345 = ((((((m01234 * a55) - (m01235 * a45)) + (m01245 * a35)) - (m01345 * a25)) + (m02345 * a15)) - (m12345 * a05));
    return m012345;
}


inline int orientationC5( double p0, double p1, double p2, double p3, double p4, double q0, double q1, double q2, double q3, double q4, double r0, double r1, double r2, double r3, double r4, double s0, double s1, double s2, double s3, double s4, double t0, double t1, double t2, double t3, double t4, double u0, double u1, double u2, double u3, double u4) {
    double det;
    double determinant_with_1_in_row_0_return_value;
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
    double m12;
    m12 = (r0 - q0);
    double m13;
    m13 = (s0 - q0);
    double m14;
    m14 = (t0 - q0);
    double m15;
    m15 = (u0 - q0);
    double m23;
    m23 = (s0 - r0);
    double m24;
    m24 = (t0 - r0);
    double m25;
    m25 = (u0 - r0);
    double m34;
    m34 = (t0 - s0);
    double m35;
    m35 = (u0 - s0);
    double m45;
    m45 = (u0 - t0);
    double m012;
    m012 = (((m01 * r1) - (m02 * q1)) + (m12 * p1));
    double m013;
    m013 = (((m01 * s1) - (m03 * q1)) + (m13 * p1));
    double m014;
    m014 = (((m01 * t1) - (m04 * q1)) + (m14 * p1));
    double m015;
    m015 = (((m01 * u1) - (m05 * q1)) + (m15 * p1));
    double m023;
    m023 = (((m02 * s1) - (m03 * r1)) + (m23 * p1));
    double m024;
    m024 = (((m02 * t1) - (m04 * r1)) + (m24 * p1));
    double m025;
    m025 = (((m02 * u1) - (m05 * r1)) + (m25 * p1));
    double m034;
    m034 = (((m03 * t1) - (m04 * s1)) + (m34 * p1));
    double m035;
    m035 = (((m03 * u1) - (m05 * s1)) + (m35 * p1));
    double m045;
    m045 = (((m04 * u1) - (m05 * t1)) + (m45 * p1));
    double m123;
    m123 = (((m12 * s1) - (m13 * r1)) + (m23 * q1));
    double m124;
    m124 = (((m12 * t1) - (m14 * r1)) + (m24 * q1));
    double m125;
    m125 = (((m12 * u1) - (m15 * r1)) + (m25 * q1));
    double m134;
    m134 = (((m13 * t1) - (m14 * s1)) + (m34 * q1));
    double m135;
    m135 = (((m13 * u1) - (m15 * s1)) + (m35 * q1));
    double m145;
    m145 = (((m14 * u1) - (m15 * t1)) + (m45 * q1));
    double m234;
    m234 = (((m23 * t1) - (m24 * s1)) + (m34 * r1));
    double m235;
    m235 = (((m23 * u1) - (m25 * s1)) + (m35 * r1));
    double m245;
    m245 = (((m24 * u1) - (m25 * t1)) + (m45 * r1));
    double m345;
    m345 = (((m34 * u1) - (m35 * t1)) + (m45 * s1));
    double m0123;
    m0123 = ((((m012 * s2) - (m013 * r2)) + (m023 * q2)) - (m123 * p2));
    double m0124;
    m0124 = ((((m012 * t2) - (m014 * r2)) + (m024 * q2)) - (m124 * p2));
    double m0125;
    m0125 = ((((m012 * u2) - (m015 * r2)) + (m025 * q2)) - (m125 * p2));
    double m0134;
    m0134 = ((((m013 * t2) - (m014 * s2)) + (m034 * q2)) - (m134 * p2));
    double m0135;
    m0135 = ((((m013 * u2) - (m015 * s2)) + (m035 * q2)) - (m135 * p2));
    double m0145;
    m0145 = ((((m014 * u2) - (m015 * t2)) + (m045 * q2)) - (m145 * p2));
    double m0234;
    m0234 = ((((m023 * t2) - (m024 * s2)) + (m034 * r2)) - (m234 * p2));
    double m0235;
    m0235 = ((((m023 * u2) - (m025 * s2)) + (m035 * r2)) - (m235 * p2));
    double m0245;
    m0245 = ((((m024 * u2) - (m025 * t2)) + (m045 * r2)) - (m245 * p2));
    double m0345;
    m0345 = ((((m034 * u2) - (m035 * t2)) + (m045 * s2)) - (m345 * p2));
    double m1234;
    m1234 = ((((m123 * t2) - (m124 * s2)) + (m134 * r2)) - (m234 * q2));
    double m1235;
    m1235 = ((((m123 * u2) - (m125 * s2)) + (m135 * r2)) - (m235 * q2));
    double m1245;
    m1245 = ((((m124 * u2) - (m125 * t2)) + (m145 * r2)) - (m245 * q2));
    double m1345;
    m1345 = ((((m134 * u2) - (m135 * t2)) + (m145 * s2)) - (m345 * q2));
    double m2345;
    m2345 = ((((m234 * u2) - (m235 * t2)) + (m245 * s2)) - (m345 * r2));
    double m01234;
    m01234 = (((((m0123 * t3) - (m0124 * s3)) + (m0134 * r3)) - (m0234 * q3)) + (m1234 * p3));
    double m01235;
    m01235 = (((((m0123 * u3) - (m0125 * s3)) + (m0135 * r3)) - (m0235 * q3)) + (m1235 * p3));
    double m01245;
    m01245 = (((((m0124 * u3) - (m0125 * t3)) + (m0145 * r3)) - (m0245 * q3)) + (m1245 * p3));
    double m01345;
    m01345 = (((((m0134 * u3) - (m0135 * t3)) + (m0145 * s3)) - (m0345 * q3)) + (m1345 * p3));
    double m02345;
    m02345 = (((((m0234 * u3) - (m0235 * t3)) + (m0245 * s3)) - (m0345 * r3)) + (m2345 * p3));
    double m12345;
    m12345 = (((((m1234 * u3) - (m1235 * t3)) + (m1245 * s3)) - (m1345 * r3)) + (m2345 * q3));
    double m012345;
    m012345 = ((((((m01234 * u4) - (m01235 * t4)) + (m01245 * s4)) - (m01345 * r4)) + (m02345 * q4)) - (m12345 * p4));
    determinant_with_1_in_row_0_return_value = m012345;
    det = determinant_with_1_in_row_0_return_value;
    int int_tmp_result;
    double eps;
    double max1 = fabs(p1);
    if( (max1 < fabs(q1)) )
    {
        max1 = fabs(q1);
    }
    if( (max1 < fabs(r1)) )
    {
        max1 = fabs(r1);
    }
    if( (max1 < fabs(s1)) )
    {
        max1 = fabs(s1);
    }
    if( (max1 < fabs(t1)) )
    {
        max1 = fabs(t1);
    }
    if( (max1 < fabs(u1)) )
    {
        max1 = fabs(u1);
    }
    double max2 = fabs(p2);
    if( (max2 < fabs(q2)) )
    {
        max2 = fabs(q2);
    }
    if( (max2 < fabs(r2)) )
    {
        max2 = fabs(r2);
    }
    if( (max2 < fabs(s2)) )
    {
        max2 = fabs(s2);
    }
    if( (max2 < fabs(t2)) )
    {
        max2 = fabs(t2);
    }
    if( (max2 < fabs(u2)) )
    {
        max2 = fabs(u2);
    }
    double max3 = fabs(p3);
    if( (max3 < fabs(q3)) )
    {
        max3 = fabs(q3);
    }
    if( (max3 < fabs(r3)) )
    {
        max3 = fabs(r3);
    }
    if( (max3 < fabs(s3)) )
    {
        max3 = fabs(s3);
    }
    if( (max3 < fabs(t3)) )
    {
        max3 = fabs(t3);
    }
    if( (max3 < fabs(u3)) )
    {
        max3 = fabs(u3);
    }
    double max4 = fabs(p4);
    if( (max4 < fabs(q4)) )
    {
        max4 = fabs(q4);
    }
    if( (max4 < fabs(r4)) )
    {
        max4 = fabs(r4);
    }
    if( (max4 < fabs(s4)) )
    {
        max4 = fabs(s4);
    }
    if( (max4 < fabs(t4)) )
    {
        max4 = fabs(t4);
    }
    if( (max4 < fabs(u4)) )
    {
        max4 = fabs(u4);
    }
    double max5 = fabs(m01);
    if( (max5 < fabs(m05)) )
    {
        max5 = fabs(m05);
    }
    if( (max5 < fabs(m04)) )
    {
        max5 = fabs(m04);
    }
    if( (max5 < fabs(m03)) )
    {
        max5 = fabs(m03);
    }
    if( (max5 < fabs(m02)) )
    {
        max5 = fabs(m02);
    }
    if( (max5 < fabs(m15)) )
    {
        max5 = fabs(m15);
    }
    if( (max5 < fabs(m14)) )
    {
        max5 = fabs(m14);
    }
    if( (max5 < fabs(m13)) )
    {
        max5 = fabs(m13);
    }
    if( (max5 < fabs(m12)) )
    {
        max5 = fabs(m12);
    }
    if( (max5 < fabs(m34)) )
    {
        max5 = fabs(m34);
    }
    if( (max5 < fabs(m25)) )
    {
        max5 = fabs(m25);
    }
    if( (max5 < fabs(m24)) )
    {
        max5 = fabs(m24);
    }
    if( (max5 < fabs(m23)) )
    {
        max5 = fabs(m23);
    }
    if( (max5 < fabs(m45)) )
    {
        max5 = fabs(m45);
    }
    if( (max5 < fabs(m35)) )
    {
        max5 = fabs(m35);
    }
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
    if( (lower_bound_1 < 8.19482853969781542511e-60) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 3.21387608851797912384e+60) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (6.02067348555779570000e-13 * ((((max5 * max1) * max2) * max3) * max4));
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



