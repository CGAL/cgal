
inline double determinant( double a00, double a01, double a02, double a03, double a04, double a10, double a11, double a12, double a13, double a14, double a20, double a21, double a22, double a23, double a24, double a30, double a31, double a32, double a33, double a34, double a40, double a41, double a42, double a43, double a44) {
    double m01;
    m01 = ((a10 * a01) - (a00 * a11));
    double m02;
    m02 = ((a20 * a01) - (a00 * a21));
    double m03;
    m03 = ((a30 * a01) - (a00 * a31));
    double m04;
    m04 = ((a40 * a01) - (a00 * a41));
    double m12;
    m12 = ((a20 * a11) - (a10 * a21));
    double m13;
    m13 = ((a30 * a11) - (a10 * a31));
    double m14;
    m14 = ((a40 * a11) - (a10 * a41));
    double m23;
    m23 = ((a30 * a21) - (a20 * a31));
    double m24;
    m24 = ((a40 * a21) - (a20 * a41));
    double m34;
    m34 = ((a40 * a31) - (a30 * a41));
    double m012;
    m012 = (((m12 * a02) - (m02 * a12)) + (m01 * a22));
    double m013;
    m013 = (((m13 * a02) - (m03 * a12)) + (m01 * a32));
    double m014;
    m014 = (((m14 * a02) - (m04 * a12)) + (m01 * a42));
    double m023;
    m023 = (((m23 * a02) - (m03 * a22)) + (m02 * a32));
    double m024;
    m024 = (((m24 * a02) - (m04 * a22)) + (m02 * a42));
    double m034;
    m034 = (((m34 * a02) - (m04 * a32)) + (m03 * a42));
    double m123;
    m123 = (((m23 * a12) - (m13 * a22)) + (m12 * a32));
    double m124;
    m124 = (((m24 * a12) - (m14 * a22)) + (m12 * a42));
    double m134;
    m134 = (((m34 * a12) - (m14 * a32)) + (m13 * a42));
    double m234;
    m234 = (((m34 * a22) - (m24 * a32)) + (m23 * a42));
    double m0123;
    m0123 = ((((m123 * a03) - (m023 * a13)) + (m013 * a23)) - (m012 * a33));
    double m0124;
    m0124 = ((((m124 * a03) - (m024 * a13)) + (m014 * a23)) - (m012 * a43));
    double m0134;
    m0134 = ((((m134 * a03) - (m034 * a13)) + (m014 * a33)) - (m013 * a43));
    double m0234;
    m0234 = ((((m234 * a03) - (m034 * a23)) + (m024 * a33)) - (m023 * a43));
    double m1234;
    m1234 = ((((m234 * a13) - (m134 * a23)) + (m124 * a33)) - (m123 * a43));
    double m01234;
    m01234 = (((((m1234 * a04) - (m0234 * a14)) + (m0134 * a24)) - (m0124 * a34)) + (m0123 * a44));
    return m01234;
}


inline int orientationC5( double p0, double p1, double p2, double p3, double p4, double q0, double q1, double q2, double q3, double q4, double r0, double r1, double r2, double r3, double r4, double s0, double s1, double s2, double s3, double s4, double t0, double t1, double t2, double t3, double t4, double u0, double u1, double u2, double u3, double u4) {
    double pq0;
    pq0 = (q0 - p0);
    double pq1;
    pq1 = (q1 - p1);
    double pq2;
    pq2 = (q2 - p2);
    double pq3;
    pq3 = (q3 - p3);
    double pq4;
    pq4 = (q4 - p4);
    double pr0;
    pr0 = (r0 - p0);
    double pr1;
    pr1 = (r1 - p1);
    double pr2;
    pr2 = (r2 - p2);
    double pr3;
    pr3 = (r3 - p3);
    double pr4;
    pr4 = (r4 - p4);
    double ps0;
    ps0 = (s0 - p0);
    double ps1;
    ps1 = (s1 - p1);
    double ps2;
    ps2 = (s2 - p2);
    double ps3;
    ps3 = (s3 - p3);
    double ps4;
    ps4 = (s4 - p4);
    double pt0;
    pt0 = (t0 - p0);
    double pt1;
    pt1 = (t1 - p1);
    double pt2;
    pt2 = (t2 - p2);
    double pt3;
    pt3 = (t3 - p3);
    double pt4;
    pt4 = (t4 - p4);
    double pu0;
    pu0 = (t0 - p0);
    double pu1;
    pu1 = (u1 - p1);
    double pu2;
    pu2 = (u2 - p2);
    double pu3;
    pu3 = (u3 - p3);
    double pu4;
    pu4 = (u4 - p4);
    double det;
    det = determinant( pq0, pq1, pq2, pq3, pq4, pr0, pr1, pr2, pr3, pr4, ps0, ps1, ps2, ps3, ps4, pt0, pt1, pt2, pt3, pt4, pu0, pu1, pu2, pu3, pu4 );
    int int_tmp_result;
    double eps;
    double max1 = fabs(pq0);
    if( (max1 < fabs(pr0)) )
    {
        max1 = fabs(pr0);
    }
    if( (max1 < fabs(ps0)) )
    {
        max1 = fabs(ps0);
    }
    if( (max1 < fabs(pt0)) )
    {
        max1 = fabs(pt0);
    }
    if( (max1 < fabs(pu0)) )
    {
        max1 = fabs(pu0);
    }
    double max2 = fabs(pq1);
    if( (max2 < fabs(pr1)) )
    {
        max2 = fabs(pr1);
    }
    if( (max2 < fabs(ps1)) )
    {
        max2 = fabs(ps1);
    }
    if( (max2 < fabs(pt1)) )
    {
        max2 = fabs(pt1);
    }
    if( (max2 < fabs(pu1)) )
    {
        max2 = fabs(pu1);
    }
    double max3 = fabs(pq2);
    if( (max3 < fabs(pr2)) )
    {
        max3 = fabs(pr2);
    }
    if( (max3 < fabs(ps2)) )
    {
        max3 = fabs(ps2);
    }
    if( (max3 < fabs(pt2)) )
    {
        max3 = fabs(pt2);
    }
    if( (max3 < fabs(pu2)) )
    {
        max3 = fabs(pu2);
    }
    double max4 = fabs(pq3);
    if( (max4 < fabs(pr3)) )
    {
        max4 = fabs(pr3);
    }
    if( (max4 < fabs(ps3)) )
    {
        max4 = fabs(ps3);
    }
    if( (max4 < fabs(pt3)) )
    {
        max4 = fabs(pt3);
    }
    if( (max4 < fabs(pu3)) )
    {
        max4 = fabs(pu3);
    }
    double max5 = fabs(pq4);
    if( (max5 < fabs(pr4)) )
    {
        max5 = fabs(pr4);
    }
    if( (max5 < fabs(ps4)) )
    {
        max5 = fabs(ps4);
    }
    if( (max5 < fabs(pt4)) )
    {
        max5 = fabs(pt4);
    }
    if( (max5 < fabs(pu4)) )
    {
        max5 = fabs(pu4);
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
        eps = (2.22889232457534740153e-13 * ((((max1 * max2) * max3) * max4) * max5));
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

